#ifndef _MY_DEFINES_H
#define _MY_DEFINES_H


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define NUM_PMTS		9
#define NUM_SAMPL		79
#define SAMPL_PITCH		((float) 1.50)

#define CAMERA_SIZE		((float) (NUM_SAMPL * SAMPL_PITCH))
#define CAMERA_MIN_POS		((float) (-CAMERA_SIZE / 2.00f))
#define CAMERA_MAX_POS		((float) (+CAMERA_SIZE / 2.00f))

#define SIZE_CONTR_GRID		6
#define CONTR_FACTOR		((float) 1.75)
#define NUM_CONTR_GRID_ITER	12

#define MX			3
#define MY			3
#define KX			10
#define KY			10


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The apparently useless terms of the form "0 * N" are used
// to force the compiler to parse unused inputs. Dude, chill
// out!  The compiler will then optimize those terms away!
#define MAP_2D(__dimx, __dimy, __x, __y)		((__y) * (__dimx) + (__x) + (0 * (__dimy)))
#define UNMAP_2D_X(__dimx, __dimy, __index)		(((__index) % (__dimx)) + (0 * (__dimy)))
#define UNMAP_2D_Y(__dimx, __dimy, __index)		(((__index) / (__dimx)) + (0 * (__dimy)))


// The apparently useless terms of the form "0 * N" are
// used to force the compiler to parse unused inputs.
// The compiler will then optimize those terms away!
#define MAP_3D(__dimx, __dimy, __dimz, __x, __y, __z)	(((__z) * (__dimy) + (__y)) * (__dimx) + (__x) + (0 * (__dimz)))
#define UNMAP_3D_X(__dimx, __dimy, __dimz, __index)	(((__index) % (__dimx)) + (0 * (__dimy) * (__dimz)))
#define UNMAP_3D_Y(__dimx, __dimy, __dimz, __index)	((((__index) / (__dimx)) % (__dimy)) + (0 * (__dimz)))
#define UNMAP_3D_Z(__dimx, __dimy, __dimz, __index)	(((__index) / ((__dimx) * (__dimy))) + (0 * (__dimz)))


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif // _MY_DEFINES_H
