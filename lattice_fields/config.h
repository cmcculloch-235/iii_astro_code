#ifndef CONFIG_H_INC
#define CONFIG_H_INC
/* Horrible way to adjust parameters... */
#define TRUE 1
#define FALSE 0

#define PARAM_SMOOTH_LEN 20.0
#define PARAM_SMOOTH TRUE
#define PARAM_PARALLEL_NL FALSE
#define PARAM_FG_DEBUG FALSE
#define PARAM_VARIANCE FALSE
#define PARAM_PSPEC_MODE TRUE

/* Transformations of the field to correlate */
#define PARAM_CORR_F1 transform_identity
#define PARAM_CORR_F2 transform_log_normal_contrast
#endif
