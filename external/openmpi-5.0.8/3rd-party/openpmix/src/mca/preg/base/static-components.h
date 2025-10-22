/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_preg_compress_component;
extern const pmix_mca_base_component_t pmix_mca_preg_native_component;
extern const pmix_mca_base_component_t pmix_mca_preg_raw_component;

const pmix_mca_base_component_t *pmix_mca_preg_base_static_components[] = {
  &pmix_mca_preg_compress_component, 
  &pmix_mca_preg_native_component, 
  &pmix_mca_preg_raw_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

