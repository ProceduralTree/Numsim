/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const mca_base_component_t mca_smsc_cma_component;

const mca_base_component_t *mca_smsc_base_static_components[] = {
  &mca_smsc_cma_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

