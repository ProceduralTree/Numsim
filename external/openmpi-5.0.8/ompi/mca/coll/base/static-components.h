/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const mca_base_component_t mca_coll_adapt_component;
extern const mca_base_component_t mca_coll_basic_component;
extern const mca_base_component_t mca_coll_han_component;
extern const mca_base_component_t mca_coll_inter_component;
extern const mca_base_component_t mca_coll_libnbc_component;
extern const mca_base_component_t mca_coll_self_component;
extern const mca_base_component_t mca_coll_sync_component;
extern const mca_base_component_t mca_coll_tuned_component;
extern const mca_base_component_t mca_coll_ftagree_component;
extern const mca_base_component_t mca_coll_monitoring_component;
extern const mca_base_component_t mca_coll_sm_component;

const mca_base_component_t *mca_coll_base_static_components[] = {
  &mca_coll_adapt_component, 
  &mca_coll_basic_component, 
  &mca_coll_han_component, 
  &mca_coll_inter_component, 
  &mca_coll_libnbc_component, 
  &mca_coll_self_component, 
  &mca_coll_sync_component, 
  &mca_coll_tuned_component, 
  &mca_coll_ftagree_component, 
  &mca_coll_monitoring_component, 
  &mca_coll_sm_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

