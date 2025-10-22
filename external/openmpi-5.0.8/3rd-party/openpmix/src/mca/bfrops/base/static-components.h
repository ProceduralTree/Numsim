/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_bfrops_v12_component;
extern const pmix_mca_base_component_t pmix_mca_bfrops_v20_component;
extern const pmix_mca_base_component_t pmix_mca_bfrops_v21_component;
extern const pmix_mca_base_component_t pmix_mca_bfrops_v3_component;
extern const pmix_mca_base_component_t pmix_mca_bfrops_v4_component;
extern const pmix_mca_base_component_t pmix_mca_bfrops_v41_component;

const pmix_mca_base_component_t *pmix_mca_bfrops_base_static_components[] = {
  &pmix_mca_bfrops_v12_component, 
  &pmix_mca_bfrops_v20_component, 
  &pmix_mca_bfrops_v21_component, 
  &pmix_mca_bfrops_v3_component, 
  &pmix_mca_bfrops_v4_component, 
  &pmix_mca_bfrops_v41_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

