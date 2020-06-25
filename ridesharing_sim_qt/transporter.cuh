# pragma once

#include "transporter.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


offer best_offer_all_buses_cuda();


__global__ offer best_offer_all_buses_cuda_kernel(
	int** stop_num_array,
	int** stop_times
	);