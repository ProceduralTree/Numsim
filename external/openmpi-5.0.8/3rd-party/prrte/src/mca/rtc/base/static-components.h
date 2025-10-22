/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_rtc_hwloc_component;

const pmix_mca_base_component_t *prte_rtc_base_static_components[] = {
  &prte_mca_rtc_hwloc_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

