#include <iostream>
#include <cstdint>
#include <string>
#include <future>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <parallel/algorithm>

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

#include <assert.h>

#include "include/kernels/warp_znorm.hpp"    // cpu_znorm call
#include "include/binary_IO.hpp"
#include "include/hpc_helpers.hpp"
#include "include/cbf_generator.hpp"

//#define DEBUG

#define DATABASE (0)
#define STREAM   (1)

#define TIMERSTART_CUDA(label)                                                 \
        cudaSetDevice(0);                                                      \
        cudaEvent_t start##label, stop##label;                                 \
        float time##label;                                                     \
        cudaEventCreate(&start##label);                                        \
        cudaEventCreate(&stop##label);                                         \
        cudaEventRecord(start##label, 0);

#define TIMERSTOP_CUDA(label)                                                  \
        cudaSetDevice(0);                                                      \
        cudaEventRecord(stop##label, 0);                                       \
        cudaEventSynchronize(stop##label);                                     \
        cudaEventElapsedTime(&time##label, start##label, stop##label);         \
        if (query_type == DATABASE)                                            \
            std::cerr << "TIMING: " << time##label << " ms " << ((num_features+1)*(num_features+1)*num_entries*num_queries)/(time##label*1e6) << " GCPUS (" << #label << ")" << std::endl; \
        if (query_type == STREAM)                                              \
            std::cerr << "TIMING: " << time##label << " ms " << ((num_features+1)*(num_features+1)*(num_entries-num_features+1)*num_queries)/(time##label*1e6) << " GCPUS (" << #label << ")" << std::endl;

typedef float value_t;                              // data type for values
typedef uint64_t index_t;                           // data type for indices
typedef uint8_t  label_t;                           // data type for label

// maximum number of features fitting into constant memory
constexpr index_t max_features = (1UL<<16)/sizeof(value_t);
__constant__ value_t cQuery[max_features];

#define LOCAL_ZNORM_STREAM   // <- this triggers the local znorm in ecg 1023 kernel

#include "include/ED.hpp"
#include "include/DTW.hpp"
using namespace FullDTW;


// add sort index block of code
using namespace std;

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
    thrust::device_vector<T> d_vec(v.begin(), v.end());
    //create index vector
    thrust::device_vector<size_t> idx(v.size());
    thrust::sequence(idx.begin(), idx.end());
    //sort
    thrust::sort_by_key(d_vec.begin(), d_vec.end(), idx.begin());
    // transform to std vector
    std::vector<size_t> idx_result(v.size());
    thrust::copy(idx.begin(), idx.end(), idx_result.begin());

    return idx_result;
}

template <typename T>
vector<size_t> sort_indexes_cpu(const vector<T> &v) {
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    // stable_sort(idx.begin(), idx.end(),
    //     [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    // https://gist.github.com/nico-zck/d166cd362c397fbc2245cf47cfbdcb44
    __gnu_parallel::sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}

template<class T>
void parallel_sort(T* data, int len, int grainsize)
{
    // Use grainsize instead of thread count so that we don't e.g.
    // spawn 4 threads just to sort 8 elements.
    if(len < grainsize)
    {
        std::sort(data, data + len, std::less<T>());
    }
    else
    {
        auto future = std::async(parallel_sort<T>, data, len/2, grainsize);

        // No need to spawn another thread just to block the calling 
        // thread which would do nothing.
        parallel_sort(data + len/2, len - len/2, grainsize);

        future.wait();

        std::inplace_merge(data, data + len/2, data + len, std::less<T>());
    }
}

int main (int argc, char * argv[]) {

    // 4 inputs plus ./main command == 5
    if (argc != 4) {
        std::cerr << "need 3 inputs" << "\n";
        return 0;
    }

    /********/
    index_t num_gpus = 1;                          // number of GPUs to be used
    index_t num_streams = 1;                       // number of streams per GPU
    index_t batch_size = 1UL << 17;                // size of a batch
    index_t buffer_size = num_streams*batch_size;  // total entries on one GPU
    index_t num_queries = 1;
    index_t num_entries;            //10757835; //15000000;
    index_t num_features;           //511; //255;

    // choose whether this is a database query (0) or stream query (>0)
    const uint8_t query_type = STREAM; //DATABASE; //STREAM;

    // configure working modes
    const bool subwarp = true; //false;    // true;                 // use subwarp kernels
    const bool enable_omp __attribute__((unused)) = false;  // enable/disable openmp in check
    const bool normalize_stream = false;                    // normalize queries in stream
    const bool lower_bound_stream = false;                  // apply lower_bounds to stream
    value_t * bsf = nullptr;

    // create systemwide pointer and initialize with infinity
    if (lower_bound_stream) {
        cudaMallocManaged(&bsf, sizeof(value_t));
        *bsf = INFINITY;
    } CUERR

    const bool init_db = true;                     // initialize DB with CBF
    const bool perform_memcpy = true;              // perform data transfers

    int counter = 0;
    int num_of_top_matches = atoi(argv[1]);  //10

    // initialize reference sequence
    value_t * data_cpu  = nullptr;
    std::string filename_reference = argv[2];
    
    std::ifstream ifile_ref(filename_reference.c_str(), std::ios::binary | std::ios::ate);
    std::streamsize size = ifile_ref.tellg();
    ifile_ref.seekg(0, std::ios::beg);
    
    cudaMallocHost(&data_cpu, size);         CUERR
    ifile_ref.read((char*) data_cpu, size);

    num_entries = size/sizeof(value_t);
    //std::cerr << "Size of Ref: " << num_entries << std::endl;


    // initialize query sequence
    value_t * starting_size = nullptr;
    cudaMallocHost(&starting_size, sizeof(value_t)*1);
    std::string filename_query = argv[3];
    std::ifstream ifile(filename_query.c_str(), std::ios::binary);

    
    while (ifile.peek() != EOF) {

    #ifdef DEBUG    
        TIMERSTART(generate_data)
    #endif

        // Load query
        value_t * query_cpu = nullptr;
        ifile.read((char*) starting_size, sizeof(value_t)*1);
        // std::cerr << "start size: " << starting_size[0] << std::endl;

        num_features = starting_size[0]-1;

        cudaMallocHost(&query_cpu, sizeof(value_t)*starting_size[0]);         CUERR
        ifile.read((char*) query_cpu, sizeof(value_t)*starting_size[0]);


        for (index_t gpu = 0; gpu < num_gpus; gpu++) {
            cudaSetDevice(gpu);
            cudaMemcpyToSymbol(cQuery, query_cpu,
                            sizeof(value_t)*num_features);
        } CUERR

    #ifdef DEBUG
        TIMERSTOP(generate_data)
    #endif

    #ifdef DEBUG
        TIMERSTART(malloc)
    #endif

        // some consistency checks
        assert(query_type == DATABASE || batch_size >= num_features);
        assert(num_features <= max_features);
        assert(num_queries == 1);

        // status
        const value_t CU = (num_entries-num_features+1)
                            *  num_features*num_features;
    #ifdef DEBUG
        std::cerr << "We are going to process "
                << CU/1000000000000.0
                << " Tera Cell Updates (TCU)"
                << std::endl;
    #endif
        const value_t DM = (num_entries-num_features+1)
                        * 2*sizeof(value_t);
    #ifdef DEBUG
        std::cerr << "We are going to stream at least "
                << DM/1073741824.0
                << " Gibi Bytes (GiB) to and from the GPU"
                << std::endl;
    #endif
        // create the streams on each GPU
        cudaStream_t streams[num_gpus][num_streams];
        for (index_t gpu = 0; gpu < num_gpus; gpu++) {
            cudaSetDevice(gpu);
            for (index_t stream = 0; stream < num_streams; stream++) {
                cudaStreamCreate(&streams[gpu][stream]);
            }
        }
        CUERR

        // initialize device memory
        value_t * dist_cpu  = nullptr,                 // distance array on the CPU
                * data_gpu[num_gpus],                  // buffers on GPUs
                * dist_gpu[num_gpus];                  // distance arrays on GPUs

        
        cudaMallocHost(&dist_cpu, sizeof(value_t)*(num_entries-num_features+1));
        for (index_t gpu = 0; gpu < num_gpus; gpu++) {
            cudaSetDevice(gpu);
            cudaMalloc(&data_gpu[gpu], sizeof(value_t)*buffer_size);
            cudaMalloc(&dist_gpu[gpu],
                    sizeof(value_t)*num_streams*(batch_size-num_features+1));
        } CUERR

    #ifdef DEBUG
        TIMERSTOP(malloc)
    #endif

    #ifdef DEBUG
        TIMERSTART_CUDA(streamed_computation)
    #endif

        for (index_t batch = 0; /*no a priori bound check possible*/ ; batch++) {

            // determine gpu and stream identifier from batch identifier
            const index_t gpu = batch % num_gpus;
            const index_t stream = (batch/num_gpus) % num_streams;
            cudaSetDevice(gpu);

            // range_size == batch_size in DB case but shortened by num_features
            // to account for overlap in the stream case
            const index_t range_size = query_type == DATABASE ? batch_size:
                                    batch_size-num_features+1;

            // slice the corresponding range from host memory
            const index_t lower = std::min(batch*range_size, num_entries);
            const index_t upper = std::min(lower+batch_size, num_entries);
            const index_t width = upper-lower;

            // if empty batch then exit
            if (width == 0)
                break;
            // if not enough points in last batch of stream then exit
            if (query_type == STREAM && width < num_features)
                break;

            // compute host and device pointers
            const index_t multiplicator = query_type == DATABASE ? num_features : 1;
            const auto data_ptr_gpu = data_gpu[gpu]+range_size*stream*multiplicator;
            const auto data_ptr_cpu = data_cpu     +range_size*batch*multiplicator;
            const auto dist_ptr_gpu = dist_gpu[gpu]+range_size*stream*num_queries;
            const auto dist_ptr_cpu = dist_cpu     +range_size*batch*num_queries;

            // toggle between width many time series of length num_features to be
            // copied in the DB case and width many data points in the stream case
            const index_t num_entries_data = query_type == DATABASE ?
                                            width*num_features :
                                            width;
            const index_t num_entries_dist = query_type == DATABASE?
                                            width :
                                            width-num_features+1;

            // reset score values on the GPU to 0
            cudaMemsetAsync(dist_ptr_gpu, 0,
                            sizeof(value_t)*num_entries_dist*num_queries,
                            streams[gpu][stream]);

            // copy the database batch to the GPU
            if (perform_memcpy)
                cudaMemcpyAsync(data_ptr_gpu, data_ptr_cpu,
                                sizeof(value_t)*num_entries_data,
                                cudaMemcpyHostToDevice,
                                streams[gpu][stream]);

            // here we call the distance function
            dist(data_ptr_gpu, dist_ptr_gpu,
                width, num_features, num_queries, subwarp,
            query_type, lower_bound_stream, bsf, streams[gpu][stream]);

            // copy distances back to CPU
            if (perform_memcpy)
                cudaMemcpyAsync(dist_ptr_cpu, dist_ptr_gpu,
                                sizeof(value_t)*num_entries_dist,
                                cudaMemcpyDeviceToHost,
                                streams[gpu][stream]);

        } CUERR

        // synchronize all streams
        for (index_t gpu = 0; gpu < num_gpus; gpu++) {
            cudaSetDevice(gpu);
            for (index_t stream = 0; stream < num_streams; stream++) {
                cudaStreamSynchronize(streams[gpu][stream]);
            }
        } CUERR

    #ifdef DEBUG
        TIMERSTOP_CUDA(streamed_computation)
    #endif

        if (lower_bound_stream)
            std::cerr << "STATUS: value stored in bsf: " << *bsf << std::endl;
        
    #ifdef DEBUG
        TIMERSTART(sort)
    #endif

        std::vector<float> v(dist_cpu, dist_cpu + (num_entries-num_features+1));
        vector<size_t> idx = sort_indexes(v);

        for ( unsigned int a = 0; a < num_of_top_matches; a = a + 1 ) {
            std::cout << counter << " " << idx[a] << " " << dist_cpu[idx[a]] << std::endl;
        }
        //for ( unsigned int a = 0; a < (num_entries-num_features+1); a = a + 1 ) {
        //    std::cout << counter << " " << a << " " << dist_cpu[a] << std::endl;
        //}
    #ifdef DEBUG
        TIMERSTOP(sort)
    #endif

    #ifdef DEBUG
        TIMERSTART(free)
    #endif
        // tear down all streams and GPU memory
        for (index_t gpu = 0; gpu < num_gpus; gpu++) {
            cudaSetDevice(gpu);
            for (index_t stream = 0; stream < num_streams; stream++)
                cudaStreamDestroy(streams[gpu][stream]);
            cudaFree(data_gpu[gpu]);
            cudaFree(dist_gpu[gpu]);
        } CUERR

        if (lower_bound_stream)
            cudaFree(bsf);                                                    CUERR

        // release the memory    
        cudaFreeHost(dist_cpu);                                               CUERR
        cudaFreeHost(query_cpu);                                              CUERR

    #ifdef DEBUG
        TIMERSTOP(free)
    #endif
        counter++;
    }

    ifile.close();

    // release memory
    cudaFreeHost(data_cpu);                                                   CUERR
}
