/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const mca_base_component_t mca_scoll_basic_component;
extern const mca_base_component_t mca_scoll_mpi_component;

const mca_base_component_t *mca_scoll_base_static_components[] = {
  &mca_scoll_basic_component, 
  &mca_scoll_mpi_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

