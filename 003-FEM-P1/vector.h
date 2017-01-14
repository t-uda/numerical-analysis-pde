
#ifndef VECTOR_H_20170114204300
#define VECTOR_H_20170114204300

typedef struct vector2d {
	double x;
	double y;
} Vector2D;

Vector2D vec2D_sub(Vector2D, Vector2D);
double vec2D_prod(Vector2D, Vector2D);
double length(Vector2D);

void check_allocated_or_abort(void *);
double * allocate_real_vector(size_t);

void vec_set(size_t, double *, double *);
void vec_add(size_t, double *, double *);
void vec_sub(size_t, double *, double *);
void vec_mul(size_t, double, double *);
void mat_vec_prod(size_t, double *, double **, double *);
double vec_prod(size_t, double *, double *);
double bilinear_form(size_t, double *, double **, double *);

void solve_CG(size_t, double *, double **, double *);

#endif

