#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <pthread.h>

/* Handle errors by printing an error message and exiting with a
  * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define ICAR_LATITUDE 224
#define ICAR_LONGITUDE 464
#define OBS_LATITUDE 222
#define OBS_LONGITUDE 462
#define ICAR_TIME_DIM 11322
#define NUM_YEARS 31
#define NUM_GRIDS 64
#define ICAR_DIR_1980 "/glade/p/ral/hap/trude/conus_icar/orig_and_removed_biasc_icar_data/merged_era/merged_era_hist_1980.nc"
#define OBS_DIR_1980  "/glade/p/ral/hap/common_data/Maurer_met_full/pr/nldas_met_update.obs.daily.pr.1980.nc"
#define NUM_THREADS 8

// create a struct to store dataset info
typedef struct dataset_info {
    int time_dim;
    int ncids[NUM_YEARS];
    int varids[NUM_YEARS];
    int year_lengths[NUM_YEARS];
    int lat;
    int lon;
} dataset_info_t;


// create a struct for thread arguments
typedef struct thread_arg {
    int cur_grid;
    int left_lon;
    int right_lon;
    int lower_lat;
    int upper_lat;
    float* p_icar_data;
    int* p_dif;
} thread_arg_t;

// creating lock to make sure that the array is updated
//pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

/*
given a start directory returns the dimensions of the whole dataset with precipitation values from 1980 to 2010
preconditions: takes in a directory (string)
postconditions: returns a struct that includes all the necessary information about the dataset such as dimensions, ncids, and varids
*/
dataset_info_t read_dataset(char* dir) {
    // find the length of the directory
    int dir_length = strlen(dir);

    // find the first half of this directory
    char* first_half = (char*) malloc((dir_length-6) * sizeof(char));
    memcpy(first_half, dir, dir_length-7);
   
    // variable for error checking
    int retval;
    
    // create a dataset_info
    dataset_info_t d_info;
    
    // create a variable that will store total time dimensions
    int time_total; 
    
    // create an iteration variable
    int year;
    
    for (year = 1980; year < 2011; ++year) {
        // crate a variable for the dataset ids and variable ids
        int ncid;
        int varid;
        
        // allocate a memory for the current directory
        char* cur_dir = (char*) malloc((dir_length+1) * sizeof(char));
        
        // copy the first half to the memory for the current directory
        strcpy(cur_dir, first_half);
        
        // concatenate the current year to the current_directory
        char str_year[5];
        snprintf(str_year, 10, "%d", year);
        strcat(cur_dir, str_year);
        
        // concatenate the netcdf extension
        strcat(cur_dir, ".nc");

        // import the first dataset
        if ((retval = nc_open(cur_dir, NC_NOWRITE, &ncid))) {
           ERR(retval);
        }

        // access the number of dimensions in the dataset
        int num_dim;
        if ((retval = nc_inq(ncid, &num_dim, NULL, NULL, NULL))) {
            ERR(retval);
        }

        // create an iteration variable
        int dim_id_num;

        // create the variables for dataArray dimensions
        int time;

        // find the length of each dimension and update the info struct
        for (dim_id_num = 0; dim_id_num < num_dim; ++dim_id_num) {
            size_t len_dim;
            char* name_dim = (char*) malloc(sizeof(char)*(NC_MAX_NAME+1));

            if ((retval = nc_inq_dim(ncid, dim_id_num, name_dim, &len_dim))) {
                ERR(retval);
            }

            if (strcmp(name_dim, "time") == 0) {
                if (year == 1980) {
                    d_info.time_dim = (int) len_dim;
                }
                else {
                    d_info.time_dim += (int) len_dim;
                }
                
                d_info.year_lengths[year-1980] = (int) len_dim;
            }
            else if ((strcmp(name_dim, "latitude") == 0) || (strcmp(name_dim, "lat") == 0)){
                d_info.lat = (int) len_dim;
            }
            else if ((strcmp(name_dim, "longitude") == 0) || (strcmp(name_dim, "lon") == 0)){
                d_info.lon = (int) len_dim;
            }

            free(name_dim);
        }

        // find the id for the precipitation variable
        if ((retval = nc_inq_varid(ncid, "pr", &varid))) {
            if ((retval = nc_inq_varid(ncid, "icar_pcp", &varid)))
               ERR(retval);
        }
        
        // add the ncids and varids to the info card
        d_info.ncids[year-1980] = ncid;
        d_info.varids[year-1980] = varid;
        
        // free the memory reserved for the current directory
        free(cur_dir);
    }
    
    // now that we are done with collecting information free the first_half of current directory
    free(first_half);
    
    return d_info;
}       


// given the info about different datasets, copies the data in netcdf to arr
// preconditions: takes in an array of (any time size, 224 latitude, 464 longitude values), an array of dataset ids
// an array of variable ids, and a time dimension which should match up with the length of datasets id array
// postconditions: does not return anything, copies the contents of netcdf files into the input array
void copy_data_icar(float arr[][ICAR_LATITUDE][ICAR_LONGITUDE], int ncids[], int varids[], int year_lengths[]) {
    // create a variable for error check in
    int retval;
    
    // create the iteration variable
    int year_num;
    int cur_year = 0;
    int time_dim = NUM_YEARS;
    
    // iterate over each year and read it into the array
    for (year_num = 0; year_num < time_dim; ++year_num) {
        // copy the data into this array
        if ((retval = nc_get_var_float(ncids[year_num], varids[year_num], &arr[cur_year][0][0])))
           ERR(retval);
        
        cur_year += year_lengths[year_num];
    }
}    

// given the info about different datasets, copies the data in netcdf to arr
// preconditions: takes in an array of (any time size, 224 latitude, 464 longitude values), an array of dataset ids
// an array of variable ids, and a time dimension which should match up with the length of datasets id array
// postconditions: does not return anything, copies the contents of netcdf files into the input array
void copy_data_obs(float arr[][OBS_LATITUDE][OBS_LONGITUDE], int ncids[], int varids[], int year_lengths[]) {
    // create a variable for error check in
    int retval;
    
    // create the iteration variable
    int year_num;
    int cur_year = 0;
    int time_dim = NUM_YEARS;
    
    // iterate over each year and read it into the array
    for (year_num = 0; year_num < time_dim; ++year_num) {
        // copy the data into this array
        if ((retval = nc_get_var_float(ncids[year_num], varids[year_num], &arr[cur_year][0][0])))
           ERR(retval);
        
        cur_year += year_lengths[year_num];
    }
}

// find the day with the minimum difference to given day
// preconditions: takes in an array of differences, the day we are comparing, and the size of the time axis
// postconditions: returns an integer, the index that allows us to find the minimum difference
int find_min_index(float mean_difs[], int day, int total_day) {
    // create a variable for minimum index
    int min_index; 
    
    // pick the starting day
    if (day != 0) {
        min_index = 0;
    }
    else {
        min_index = 1;
    }
    
    // create an iteration variable
    int d;
    // iterate over all the mean differences
    for (d = 0; d < total_day; d++) {
        // if we are not on the day being compared and the current min difs is less than the prev min dif
        if ((day != d) && (mean_difs[d] < mean_difs[min_index])) {
            // update min index
            min_index = d;
        }
    }
    // return the minimum index
    return min_index;
}

// store all the days with closest precipitation rates in an array
// preconditions: takes in full icar data, a closest index array to update, and total time dimension
// postconditions: does not return anything, updates the integer array named dif_nw with indices
void find_difs(float* icar_data, int* dif_nw, int time_dim, int left_lon_icar, int right_lon_icar, int cur_grid, int lower_lat_icar, int upper_lat_icar) {
    // set the iteration variables 
    int cur_left_lon;
    int cur_lower_lat;
    
    for (cur_left_lon = left_lon_icar; cur_left_lon < right_lon_icar; cur_left_lon += 10) {
        for (cur_lower_lat = lower_lat_icar; cur_lower_lat < upper_lat_icar; cur_lower_lat += 10) {
            // find the other edges of this grid
            int cur_right_lon = cur_left_lon + 10;
            int cur_upper_lat = cur_lower_lat + 10;
            
            // inner loop iteration variables
            int lat;
            int lon;
            int day;
            int comp_day;
            // iterate over time
            for (day = 0; day < time_dim; ++day){
                float mean_difs_for_this_grid[time_dim];
                for (comp_day = 0; comp_day < time_dim; ++comp_day)  {
                    // store the differences for this grid and day in an array
                    float sum_difs = 0;
                    for (lat = cur_lower_lat; lat < cur_upper_lat; lat++) {
                        for (lon = cur_left_lon; lon < cur_right_lon; lon++) {
                            
                            float day_val = *(icar_data + day * (ICAR_LATITUDE * ICAR_LONGITUDE) + lat * ICAR_LONGITUDE + lon);
                            float comp_val = *(icar_data + comp_day * (ICAR_LATITUDE * ICAR_LONGITUDE) + lat * ICAR_LONGITUDE + lon);
                            // calculate the mean square difference
                            float dif = (day_val - comp_val) * (day_val - comp_val);
                            // add it to the sum
                            sum_difs += dif;
                        } 
                    }
                    // take the mean difference of this grid
                    mean_difs_for_this_grid[comp_day] = sum_difs / 100;   
                }
                // find the corresponding day that is closest to this day
                int min_index = find_min_index(mean_difs_for_this_grid, day, comp_day);
                
                int* cur_pointer = dif_nw + (cur_grid * ICAR_TIME_DIM + day);
                                    
                //pthread_mutex_lock(&lock);
                *cur_pointer = min_index;
                //pthread_mutex_unlock(&lock);

            }
            
            // when we are done, print the grid
            printf("%d\n", cur_grid);
            cur_grid += 1;
        }
    }
}

// finds all the corresponding observations in the northwest region and updates the float array named cor_obs
// preconditions: takes in the corresponding observation array to update, corresponding index array, full observation dataset,
// total number of grids, the total days we have
// postconditions: does not return anything, updates the cor_obs array and bias corrects the north west region in the full_icar dataset
void gen_corresponding_obs(float cor_obs[][80][80], int dif_nw [][ICAR_TIME_DIM], float full_obs[][OBS_LATITUDE][OBS_LONGITUDE], int num_grids, int time_dim) {
    // set the indices for the observation data to only consider the north_west region
    int left_lon_obs = 0;
    int right_lon_obs = 78;
    
    int lower_lat_obs = 119;
    int upper_lat_obs = 199;
    
    // create an iteration variable
    int g;
    int time;
    
    // iterate over each grid an time
    for (g = 0; g < num_grids; ++g) {
        // access the coordinates of this grid in cor_obs
        int low_lat_icar = 10 * (g % 8);
        int up_lat_icar = low_lat_icar + 10;
            
        int left_lon_icar = 10 * (g / 8);
        int right_lon_icar = left_lon_icar + 10;
            
        // iteration variables
        int lat;
        int lon;
        
        // iterate over time;
        for (time = 0; time < time_dim; ++time) {
            int time_index = dif_nw[g][time];
            
            // iterate over the coordinates
            for (lat = low_lat_icar; lat < up_lat_icar; lat++) {
                for (lon = left_lon_icar; lon < right_lon_icar; lon++) {
                    // pad the zeros
                    if (lon < 2) {
                        cor_obs[time][lat][lon] = 0;
                    }
                    else {
                        // assign the corresponding observation to this coordinate
                        cor_obs[time][lat][lon] = full_obs[time_index+1][119 + lat][left_lon_obs + (lon-2)];
                    }
                }
            } 
        }
        // print the grid number
        printf("%d\n", g);
    }
}

// a function to call find_difs using threads
void* thread_func(void* p) {
  thread_arg_t* arg = (thread_arg_t*) p;
  // find the differences in this section
  find_difs(arg->p_icar_data, arg->p_dif, ICAR_TIME_DIM, arg->left_lon, arg->right_lon, arg->cur_grid, arg->lower_lat, arg->upper_lat);
  return NULL;
}

int main(int argc, char** argv) {
    // Allocate a space for the dataset id and variable id
    int* ncid_try = (int*) malloc(sizeof(int));
    int* varid_try = (int*) malloc(sizeof(int));
    int retval;
        
    // Access the dimensions and ids of the 1980 icar dataset
    dataset_info_t icar_info = read_dataset(ICAR_DIR_1980); 
    
    printf("Icar Dimensions: %d, %d, %d\n", icar_info.time_dim, icar_info.lat, icar_info.lon);
    
    float full_icar[icar_info.time_dim][icar_info.lat][icar_info.lon];

    // copy the icar precipitation rate to this three dimensional array
    copy_data_icar(full_icar, icar_info.ncids, icar_info.varids, icar_info.year_lengths);
    
    float* p_icar = &full_icar[0][0][0];
    
    // create an array to store the days with closest precipiation rates
    int dif_nw[NUM_GRIDS][ICAR_TIME_DIM];
    
    int* p_dif = &dif_nw[0][0];
 
    // Access the dimensions and ids of the 1980 icar dataset
    dataset_info_t obs_info = read_dataset(OBS_DIR_1980); 
    
    printf("Observation Dimensions: %d, %d, %d\n", obs_info.time_dim, obs_info.lat, obs_info.lon);
    
    // reserve a space in the memory for the new dataset
    float full_obs[obs_info.time_dim][obs_info.lat][obs_info.lon];
    
    // copy the icar precipitation rate to this three dimensional array
    copy_data_obs(full_obs, obs_info.ncids, obs_info.varids, obs_info.year_lengths);
        
    // find the days with closest precipitation rate for each day and grid
    
    // allocate a memory for threads
    pthread_t threads[NUM_THREADS];
    
    // allocate a memory for thread arguments
    thread_arg_t args[NUM_THREADS];
    int i;
    
    for (i = 0; i < NUM_THREADS; ++i) {
        args[i].cur_grid = i * (NUM_GRIDS / NUM_THREADS);
        args[i].left_lon = ((i * 8) / NUM_THREADS) * 10;
        args[i].right_lon = args[i].left_lon + 10;
        args[i].lower_lat = 120 + 80 * (i % (NUM_THREADS / 8)) * (8.0 / (float) NUM_THREADS);
        args[i].upper_lat = args[i].lower_lat + 80 * (8.0 / (float) NUM_THREADS);
        args[i].p_icar_data = p_icar;
        args[i].p_dif = p_dif;
        
        printf("THREAD: %d, cur_grid: %d, left_lon: %d, right_lon: %d, lower_lat: %d, upper_lat: %d\n", i, args[i].cur_grid, args[i].left_lon, args[i].right_lon, args[i].lower_lat, args[i].upper_lat);
        
        if (pthread_create(&threads[i], NULL, thread_func, &args[i])) {
            perror("pthread_create failed");
            exit(2);
        }
    }
    
    // run the threads
    for (i = 0; i < NUM_THREADS; ++i) {
        if (pthread_join(threads[i], NULL)) {
          perror("pthread_join failed");
          exit(2);
        }
    }
    
    printf("-----------------------------------------------\n");
    
    // create an array to store the bias corrected output
    float nw_cor_obs[icar_info.time_dim][80][80];
    
    // find the corresponding observations and assign it to the bias corrected dataset
    gen_corresponding_obs(nw_cor_obs, dif_nw, full_obs, 64, icar_info.time_dim);
    
    // print some numbers in this dataset
    printf("%f\n", nw_cor_obs[0][40][20]);
    printf("%f\n", nw_cor_obs[2][62][15]);
    printf("%f\n", nw_cor_obs[0][50][20]);
    printf("%f\n", nw_cor_obs[0][62][15]);
    printf("%f\n", nw_cor_obs[0][32][33]);
    
    //create the netcdf file
    // Allocate space for netCDF dimension ids
    int timeId, latId, lonId;
  
    // Allocate space for the netcdf file id
    int ncid;
  
    // Allocate space for the data variable ids
    int pcp_id;
    
    // array to store the dimension ids
    int dimids[3];
    
    // create a netcdf file named "updated_C_data.nc"
    if((retval = nc_create("updated_C_data.nc", NC_NETCDF4, &ncid))) {
        ERR(retval);
    }
    
    // define the dimensions (time, latitude, and longitude)
    if((retval = nc_def_dim(ncid, "time", icar_info.time_dim, &timeId))) { 
        ERR(retval); 
    }
    
    if((retval = nc_def_dim(ncid, "latitude", 80, &latId))) { 
        ERR(retval); 
    }
     
    if((retval = nc_def_dim(ncid, "longitude", 80, &lonId))) { 
        ERR(retval); 
    }
    
    dimids[0] = timeId;
    dimids[1] = latId;
    dimids[2] = lonId;
    
    // Add the variable
    if((retval = nc_def_var(ncid, "pcp", NC_FLOAT, 3, dimids, &pcp_id))) { 
        ERR(retval); 
    }
    
    // End "Metadata" mode
    if((retval = nc_enddef(ncid))) { 
        ERR(retval); 
    }
    
    // Add the bias-corrected data to this netcdf file
    if((retval = nc_put_var(ncid, pcp_id, &nw_cor_obs[0][0][0]))) { 
        ERR(retval); 
    }
    
    // Close the net_cdf file for bias-corrected-data
    if((retval = nc_close(ncid))) { 
        ERR(retval); 
    }
    return 0;   
}

// add a make file
// limit the memory transfer
// make sure that memory accesses are successive
// read about efficient memory allocations