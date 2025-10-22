/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_psquash_flex128_component;
extern const pmix_mca_base_component_t pmix_mca_psquash_native_component;

const pmix_mca_base_component_t *pmix_mca_psquash_base_static_components[] = {
  &pmix_mca_psquash_flex128_component, 
  &pmix_mca_psquash_native_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

