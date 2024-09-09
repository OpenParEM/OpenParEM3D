#mpirun -np 12 OpenParEM3D filter_1st_order.proj >& filter_1st_order.log
#mpirun -np 12 OpenParEM3D filter_2nd_order.proj >& filter_2nd_order.log
mpirun -np 12 OpenParEM3D filter.proj >& filter.log

