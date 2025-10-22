/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const mca_base_component_t mca_hook_comm_method_component;

const mca_base_component_t *mca_hook_base_static_components[] = {
  &mca_hook_comm_method_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

