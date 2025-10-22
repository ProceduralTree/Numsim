/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_prteinstalldirs_env_component;
extern const pmix_mca_base_component_t prte_mca_prteinstalldirs_config_component;

const pmix_mca_base_component_t *prte_prteinstalldirs_base_static_components[] = {
  &prte_mca_prteinstalldirs_env_component, 
  &prte_mca_prteinstalldirs_config_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

