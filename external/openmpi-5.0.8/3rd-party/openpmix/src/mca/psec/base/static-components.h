/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_psec_native_component;
extern const pmix_mca_base_component_t pmix_mca_psec_none_component;

const pmix_mca_base_component_t *pmix_mca_psec_base_static_components[] = {
  &pmix_mca_psec_native_component, 
  &pmix_mca_psec_none_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

