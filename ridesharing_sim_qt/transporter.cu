#include "transporter.cuh"


offer best_offer_all_buses_cuda()
{

}


__global__ offer best_offer_all_buses_cuda_kernel()
{
	// if thread is pickup-thread:
	// compute pickup_is_possible array
	// AND
	// pickup_insertion times
	// AND
	// occupancy-at-pickup array

	__syncThreads();

	// if thread is dropoff-thread:
	// compute dropoff_is_possible array
	// AND
	// dropoff_insertion times
	// AND
	// occupancy-at-pickup array

	__syncThreads();


}