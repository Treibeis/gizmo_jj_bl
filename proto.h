#ifndef ALLVARS_H
#include "allvars.h"
#endif
#include "gravity/forcetree.h"
#include "domain.h"
#ifdef COOLING
#include "cooling/cooling.h"
#endif
#ifdef ADJ_BOX_POWERSPEC
#include "power_spec/adj_box_powerspec_proto.h"
#endif
#ifdef BLACK_HOLES
#include "./galaxy_sf/blackholes/blackhole.h"
#endif

/* declarations of functions throughout the code */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part (adding/removing routines as necessary) 
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef HAVE_HDF5
#include <hdf5.h>
void write_header_attributes_in_hdf5(hid_t handle);
void read_header_attributes_in_hdf5(char *fname);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_units_attributes_in_hdf5(hid_t handle);
void write_constants_attributes_in_hdf5(hid_t handle);
#endif
void report_VmRSS(void);
void output_compile_time_options(void);

void report_pinning(void);
void pin_to_core_set(void);
void get_core_set(void);

void init_turb(void);
void set_turb_ampl(void);
void add_turb_accel(void);
void reset_turb_temp(void);
void log_turb_temp(void);

void mpi_report_comittable_memory(long long BaseMem);
long long report_comittable_memory(long long *MemTotal,
                                   long long *Committed_AS,
                                   long long *SwapTotal,
                                   long long *SwapFree);

void merge_and_split_particles(void);
int does_particle_need_to_be_merged(MyIDType i);
int does_particle_need_to_be_split(MyIDType i);
double ref_mass_factor(MyIDType i);
void merge_particles_ij(MyIDType i, MyIDType j);
void split_particle_i(MyIDType i, int n_particles_split, MyIDType i_nearest, double r2_nearest);

void do_first_halfstep_kick(void);
void do_second_halfstep_kick(void);
void find_timesteps(void);
#ifdef GALSF
void compute_stellar_feedback(void);
#endif
void compute_hydro_densities_and_forces(void);
void compute_grav_accelerations(void);
void calc_memory_checksum(void *base, size_t bytes);

void get_disk_forces(double RR, double zz, double *f_R, double *f_z);
double get_disk_mass(double time);
void growing_disk_init(void);

double get_turb_pot(double x, double y, double z);

void   sub_turb_move_perturbers(double t0, double t1);
void   sub_turb_add_forces(void);
void   sub_turb_read_table(void);
void   sub_turb_parent_halo_accel(double dx, double dy, double dz, double *acc);
double sub_turb_enclosed_mass(double r, double msub, double vmax, double radvmax, double c);


int powerspec_turb_treefind(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local);
int powerspec_turb_find_nearest_evaluate(int target, int mode, int *nexport, int *nsend_local);
void powerspec_turb_calc_dispersion(void);
double powerspec_turb_obtain_fields(void);
void powerspec_turb_save(char *fname, double *disp);
void powerspec_turb_collect(void);
void powerspec_turb(int filenr);


void set_cosmo_factors_for_current_time(void);
void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

void set_non_standard_physics_for_current_time(void);
void calculate_non_standard_physics(void);
void compute_statistics(void);
void execute_resubmit_command(void);
void make_list_of_active_particles(void);
void output_extra_log_messages(void);


static inline double WRAP_POSITION_UNIFORM_BOX(double x)
{
    while(x >= All.BoxSize) {x -= All.BoxSize;}
    while(x < 0) {x += All.BoxSize;}
    return x;
}

static inline double DMAX(double a, double b) { return (a > b) ? a : b; }
static inline double DMIN(double a, double b) { return (a < b) ? a : b; }
static inline int IMAX(int a, int b) { return (a > b) ? a : b; } 
static inline int IMIN(int a, int b) { return (a < b) ? a : b; }
static inline double MINMOD(double a, double b) {return (a>0) ? ((b<0) ? 0 : DMIN(a,b)) : ((b>=0) ? 0 : DMAX(a,b));}
/* special version of MINMOD below: a is always the "preferred" choice, b the stability-required one. here we allow overshoot, just not opposite signage */
static inline double MINMOD_G(double a, double b) {return a;}


#ifdef SHEARING_BOX
void calc_shearing_box_pos_offset(void);
#endif


void do_distortion_tensor_kick(int i, double dt_gravkick);
void set_predicted_sph_quantities_for_extra_physics(int i);
void do_sph_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

void check_particle_for_temperature_minimum(int i);

double get_pressure(int i);

void read_fof(int num);
int fof_compare_ID_list_ID(const void *a, const void *b);

void myfree_msg(void *p, char *msg);
void kspace_neutrinos_init(void);

#ifdef TWOPOINT_FUNCTION_COMPUTATION_ENABLED
void twopoint(void);
void twopoint_save(void);
int twopoint_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat rsearch, int target, int *startnode, int mode, int *nexport, int *nsend_local);
int twopoint_count_local(int target, int mode, int *nexport, int *nsend_local);
#endif

void powerspec(int flag, int *typeflag);
double PowerSpec_Efstathiou(double k);
void powerspec_save(void);
void foldonitself(int *typelist);
void dump_potential(void);

int snIaheating_evaluate(int target, int mode, int *nexport, int *nSend_local);
void snIa_heating(void);
void voronoi_setup_exchange(void);

double get_neutrino_powerspec(double k, double ascale);
double get_powerspec(double k, double ascale);
void init_transfer_functions(void);


int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                       int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int MPI_Sizelimited_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
			     int dest, int sendtag, void *recvbuf, int recvcount,
			     MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical);
void sort_based_on_field(void *data, int field_offset, int n_items, int item_size, void **data2ptr);
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size);

void parallel_sort_special_P_GrNr_ID(void);
void calculate_power_spectra(int num, long long *ntot_type_all);

int pmforce_is_particle_high_res(int type, MyDouble *pos);

void compare_partitions(void);
void assign_unique_ids(void);
int permut_data_compare(const void *a, const void *b);
void  generate_permutation_in_active_list(void);
void get_particle_numbers(char *fname, int num_files);

void conduction(void);
void conduction_matrix_multiply(double *in, double *out);
double conduction_vector_multiply(double *a, double *b);
int conduction_evaluate(int target, int mode, double *in, double *out, double *sum,
			int *nexport, int *nsend_local);


void fof_get_group_center(double *cm, int gr);
void fof_get_group_velocity(double *cmvel, int gr);
int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);
void fof_compute_group_properties(int gr, int start, int len);

void parallel_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
void parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *), MPI_Comm comm);
int compare_IDs(const void *a, const void *b);
void test_id_uniqueness(void);

void check_particles_info(const char *func, const char *file, int linenr);

int io_compare_P_ID(const void *a, const void *b);
int io_compare_P_GrNr_SubNr(const void *a, const void *b);


void drift_particle(int i, integertime time1);
int ShouldWeDoDynamicUpdate(void);

void put_symbol(double t0, double t1, char c);
void write_cpu_log(void);

int get_timestep_bin(integertime ti_step);

const char* svn_version(void);

void find_particles_and_save_them(int num);
void lineofsight_output(void);
void sum_over_processors_and_normalize(void);
void absorb_along_lines_of_sight(void);
void output_lines_of_sight(int num);
integertime find_next_lineofsighttime(integertime time0);
integertime find_next_gridoutputtime(integertime ti_curr);
void add_along_lines_of_sight(void);
void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent, int mode);


void x86_fix(void) ;

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);

void mymalloc_init(void);
void dump_memory_table(void);
void report_detailed_memory_usage_of_largest_task(size_t *OldHighMarkBytes, const char *label, const char *func, const char *file, int line);

double get_shear_viscosity(int i);

void kinetic_feedback_mhm(void);
int kin_compare_key(const void *a, const void *b);
void kinetic_evaluate(int target, int mode);

void bubble(void);
void multi_bubbles(void);
void bh_bubble(double bh_dmass, MyFloat center[3], MyIDType BH_id);
void find_CM_of_biggest_group(void);
int compare_length_values(const void *a, const void *b);
double rho_dot(double z, void *params);
double bhgrowth(double z1, double z2);


int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);

double INLINE_FUNC Get_Particle_Size(int i);
double INLINE_FUNC Particle_density_for_energy_i(int i);
double INLINE_FUNC Get_Particle_Expected_Area(double h);
#ifdef COSMIC_RAYS
double INLINE_FUNC Get_Particle_CosmicRayPressure(int i);
double Get_CosmicRayGradientLength(int i);
double Get_CosmicRayStreamingVelocity(int i);
#endif
double INLINE_FUNC Particle_effective_soundspeed_i(int i);
#ifdef MAGNETIC
double INLINE_FUNC Get_Particle_BField(int i_particle_id, int k_vector_component);
double Get_DtB_FaceArea_Limiter(int i);
#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Particle_PhiField(int i_particle_id);
double INLINE_FUNC Get_Particle_PhiField_DampingTimeInv(int i_particle_id);
#endif
#endif

double INLINE_FUNC hubble_function(double a);
#ifdef DARKENERGY
double DarkEnergy_a(double);
double DarkEnergy_t(double);
#ifdef TIMEDEPDE
void fwa_init(void);
double INLINE_FUNC fwa(double);
double INLINE_FUNC get_wa(double);
#ifdef TIMEDEPGRAV
double INLINE_FUNC dHfak(double a);
double INLINE_FUNC dGfak(double a);
#endif
#endif
#endif

#ifdef EXTERNALHUBBLE
double INLINE_FUNC hubble_function_external(double a);
#endif

void blackhole_accretion(void);
int blackhole_evaluate(int target, int mode, int *nexport, int *nsend_local);
int blackhole_evaluate_swallow(int target, int mode, int *nexport, int *nsend_local);

int  blackhole_compare_key(const void *a, const void *b);


int ngb_treefind_fof_nearest(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local, int MyFOF_PRIMARY_LINK_TYPES);

#if defined(POPIII) && defined(GALSF) && defined(METALS)
double Legacy_PopIII(int i);
void InitLegacy();
double GetMass(int i);
#endif

#if defined(POPII_P2L) && defined(GALSF) && defined(METALS)
void P2L(int i, int mode);
#endif

#ifdef FOF_SCALEDEPENDENT
double get_next_FoF_time(double haloDynamicalTime);
#endif
#ifdef FOFRAD
int fof_rad(int snapnum);
#endif
#ifdef HOTGAS_QUENCH
void hotgas_quenching(void);
#endif
#ifdef GALSF_INSTANTANEOUS_METALS
double get_yield(double Zmet, int k, double age);
#if defined(GALSF_TYPEIA) || defined(GALSF_AGBFEEDBACK)
void Delayed_Feedback(void);
#endif
#endif
void fof_fof(int num);
void fof_import_ghosts(void);
void fof_course_binning(void);
void fof_find_groups(void);
void fof_check_cell(int p, int i, int j, int k);
void fof_find_minids(void);
int fof_link_accross(void);
void fof_exchange_id_lists(void);
int fof_grid_compare(const void *a, const void *b);
void fof_compile_catalogue(void);
void fof_save_groups(int num);
void fof_save_local_catalogue(int num);
void fof_find_nearest_dmparticle(void);
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *nexport, int *nsend_local);

int fof_compare_key(const void *a, const void *b);
void fof_link_special(void);
void fof_link_specialpair(int p, int s);
void fof_make_black_holes(void);

int io_compare_P_GrNr_ID(const void *a, const void *b);

void write_file(char *fname, int readTask, int lastTask);

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last);

int get_values_per_blockelement(enum iofields blocknr);

int get_datatype_in_block(enum iofields blocknr);
void get_dataset_name(enum iofields blocknr, char *buf);


int blockpresent(enum iofields blocknr);
void fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);

int get_particles_in_block(enum iofields blocknr, int *typelist);

int get_bytes_per_blockelement(enum iofields blocknr, int mode);

void read_file(char *fname, int readTask, int lastTask);

void get_Tab_IO_Label(enum iofields blocknr, char *label);


void long_range_init_regionsize(void);

int find_files(char *fname);

int metals_compare_key(const void *a, const void *b);
void enrichment_evaluate(int target, int mode);

int hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *hydro_evaluate_primary(void *p);
void *hydro_evaluate_secondary(void *p);

void pm_init_nonperiodic_allocate(void);

void  pm_init_nonperiodic_free(void);

double get_random_number(MyIDType id);
void set_random_numbers(void);

int grav_tree_compare_key(const void *a, const void *b);
int dens_compare_key(const void *a, const void *b);
int hydro_compare_key(const void *a, const void *b);

int data_index_compare(const void *a, const void *b);
int peano_compare_key(const void *a, const void *b);

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_domain(void *b, size_t n, size_t s);
void mysort_idlist(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
void mysort_peano(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));

void check_wind_creation(void);
void treat_outflowing_particles(void);
void set_injection_accel(void);


int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *density_evaluate_primary(void *p);
void *density_evaluate_secondary(void *p);
int density_isactive(int n);

size_t sizemax(size_t a, size_t b);


void reconstruct_timebins(void);

void init_peano_map(void);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
peanokey peano_and_morton_key(int x, int y, int z, int bits, peanokey *morton);
peanokey morton_key(int x, int y, int z, int bits);

void catch_abort(int sig);
void catch_fatal(int sig);
void terminate_processes(void);
void enable_core_dumps_and_fpu_exceptions(void);
void write_pid_file(void);

#ifdef PAUSE_RUN_TO_ATTACH_DEBUGGER
void pause_run_to_attach_debugger();
#endif

void pm_init_periodic_allocate(void);

void pm_init_periodic_free(void);

void move_particles(integertime time1);


void find_next_sync_point_and_drift(void);
void find_dt_displacement_constraint(double hfac);
#ifdef WAKEUP
void process_wake_ups(void);
#endif

void set_units_sfr(void);

void gravity_forcetest(void);

void allocate_commbuffers(void);
void allocate_memory(void);
void begrun(void);
void check_omega(void);
void close_outputfiles(void);
void compute_accelerations(void);
void compute_global_quantities_of_system(void);
void compute_potential(void);
void construct_timetree(void);
void cooling_and_starformation(void);
#ifdef SINKS
void do_sinks(void);
#endif

#if defined(TURB_DRIVING)
void do_turb_driving_step_first_half(void);
void do_turb_driving_step_second_half(void);
#endif

double evaluate_NH_from_GradRho(MyFloat gradrho[3], double hsml, double rho, double numngb_ndim, double include_h);
double evaluate_dL_from_GradRho(MyFloat gradrho[3], double hsml, double rho, double numngb_ndim, double include_h);

#ifdef GALSF
double evaluate_stellar_age_Gyr(double stellar_tform);
double evaluate_l_over_m_ssp(double stellar_age_in_gyr);
double calculate_relative_light_to_mass_ratio_from_imf(MyIDType i);
#endif
#if defined(GALSF_FB_RPWIND_LOCAL) && defined(GALSF_FB_RPWIND_FROMSTARS)
int ngb_treefind_newstars(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local);
#endif

#ifdef GRAIN_FLUID
void apply_grain_dragforce(void);
#ifdef GRAIN_COLLISIONS
void grain_collisions(void);
void grain_density(void);
int grain_density_evaluate(int target, int mode, int *nexport, int *nsend_local);
int grain_density_isactive(int n);
int grain_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local);
#endif
#endif

#ifdef GALSF_FB_GASRETURN
void stochastic_gas_return_singledomain(void);
#endif

#ifdef GALSF_FB_HII_HEATING
void HII_heating_singledomain(void);
#ifdef GALSF_FB_HII_HEATING_USEMULTIDOMAINSHARE
void HII_heating_withMPIcomm(void);
int HIIheating_RHIIest(int target);
int HIIheating_evaluate(int target, int mode, int *nexport, int *nsend_local);
#endif
#endif

#ifdef GALSF_FB_LOCAL_UV_HEATING
void selfshield_local_incident_uv_flux(void);
#endif

#if defined(GALSF_FB_SNE_HEATING) || defined(GALSF_FB_GASRETURN)
void snII_heating_singledomain(void);
void determine_where_SNe_occur(void);
void mechanical_fb_calc(int feedback_type);
int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int feedback_type);
void *addFB_evaluate_primary(void *p, int feedback_type);
void *addFB_evaluate_secondary(void *p, int feedback_type);
#ifdef GALSF_FB_SNE_HEATING_USEMULTIDOMAINSHARE
void snII_heating_withMPIcomm(void);
int snIIheating_evaluate(int target, int mode, int *nexport, int *nsend_local);
#endif
#ifdef GALSF_GASOLINE_RADHEATING
void luminosity_heating_gasoline(void);
#endif
#endif

#ifdef COOL_METAL_LINES_BY_SPECIES
/*double GetMetalLambda(double, double);*/
double getSpCoolTableVal(long i,long j,long k,long tblK);
double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z);
double GetLambdaSpecies(int k_species, double nHcgs, double logT, double fHe);
void LoadMultiSpeciesTables(void);
void ReadMultiSpeciesTables(int iT);
char *GetMultiSpeciesFilename(int i, int hk);
#endif

#if defined(BH_PHOTONMOMENTUM) || defined(BH_BAL_WINDS)
double bh_angleweight(double bh_lum_input, MyFloat bh_angle[3], double hR, double dx, double dy, double dz, int mode);
double bh_angleweight_localcoupling(int j, double hR, double theta);
#endif

#if defined(GALSF_FB_RPWIND_DO_IN_SFCALC) || defined(GALSF_SUBGRID_WINDS) || defined(GALSF_SUBGRID_VARIABLEVELOCITY)
void assign_wind_kick_from_sf_routine(int i, double sm, double dtime, double* pvtau_return);
#endif
#ifdef GALSF_SUBGRID_VBOOST
double wind_vel_boost(int i,double v_circ,double v);
#endif

#if defined(GALSF_FB_RPWIND_FROMSTARS) && !defined(GALSF_FB_RPWIND_DO_IN_SFCALC)
void radiation_pressure_winds_consolidated(void);
#endif

#if defined(BLACK_HOLES)
int blackhole_evaluate_PREPASS(int target, int mode, int *nexport, int *nSend_local);
#endif

#ifdef  GALSF_SUBGRID_DMDISPERSION
int disp_gravity_kernel_shared_check(short int particle_type_primary, short int particle_type_secondary);
void disp_setup_smoothinglengths(void);
int dm_disp_ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
                                       int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                       int *ngblist, int type_of_searching_particle);
void disp_density(void);
int disp_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *disp_density_evaluate_primary(void *p);
void *disp_density_evaluate_secondary(void *p);
int disp_density_isactive(MyIDType i);
#endif

#ifdef PM_HIRES_REGION_CLIPDM
int ngb_treefind_variable_threads_nongas(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
                                         int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                         int *ngblist);
#endif

void cooling_only(void);
void count_hot_phase(void);
void delete_node(int i);
void density(void);
void density_decouple(void);
void determine_interior(void);
int dissolvegas(void);
void do_box_wrapping(void);
double enclosed_mass(double R);
void endrun(int);
void energy_statistics(void);
void ensure_neighbours(void);

void output_log_messages(void);
void ewald_corr(double dx, double dy, double dz, double *fper);

void ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void ewald_force_ni(int iii, int jjj, int kkk, double x[3], double force[3]);

void ewald_init(void);
double ewald_psi(double x[3]);
double ewald_pot_corr(double dx, double dy, double dz);
int find_ancestor(int i);
integertime find_next_outputtime(integertime time);
void find_next_time(void);
integertime find_next_time_walk(int node);
void free_memory(void);
void advance_and_find_timesteps(void);
integertime get_timestep(int p, double *a, int flag);

void determine_PMinterior(void);
void gravity_tree(void);
void hydro_force(void);
void init(void);
void do_the_cooling_for_particle(int i);
double get_starformation_rate(int i);
int determine_sf_flag(int i);
void update_internalenergy_for_galsf_effective_eos(int i, double tcool, double tsfr, double x, double rateOfSF);
void init_clouds(void);
int ionized(int i);
void integrate_sfr(void);
void insert_node(int i);
int mark_targets(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void mpi_printf(const char *fmt, ...);
void open_outputfiles(void);
void write_outputfiles_header(void);
void peano_hilbert_order(void);
double pot_integrand(double xx);
void predict(double time);
void predict_collisionless_only(double time);
void predict_sph_particles(double time);
void prepare_decouple(void);
void read_ic(char *fname);
int read_outputlist(char *fname);
void read_parameter_file(char *fname);
void rearrange_particle_sequence(void);
void reorder_gas(void);
void reorder_particles(void);
void restart(int mod);
void run(void);
void savepositions(int num);
void savepositions_ioformat1(int num);
double my_second(void);
void set_softenings(void);
void set_sph_kernel(void);
void set_units(void);
void setup_smoothinglengths(void);

void minimum_large_ints(int n, long long *src, long long *res);
void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);

void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
void veldisp_ensure_neighbours(int mode);


double get_gravkick_factor(integertime time0, integertime time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double growthfactor_integ(double a, void *param);
double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(integertime time0, integertime time1);
double measure_time(void);
double report_time(void);

/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */

double pow(double, double);


void long_range_init(void);
void long_range_force(void);
void pm_init_periodic(void);
void pmforce_periodic(int mode, int *typelist);
void pm_init_regionsize(void);
void pm_init_nonperiodic(void);
int pmforce_nonperiodic(int grnr);

int pmpotential_nonperiodic(int grnr);
void pmpotential_periodic(void);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

double enclosed_mass(double R);
void pm_setup_nonperiodic_kernel(void);


#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
int rt_get_source_luminosity(MyIDType i, double sigma_0, double *lum);
double rt_kappa(MyIDType j, int k_freq);
double rt_absorption_rate(MyIDType i, int k_freq);
double rt_diffusion_coefficient(MyIDType i, int k_freq);
void rt_eddington_update_calculation(MyIDType j);
void rt_update_driftkick(MyIDType i, double dt_entr, int mode);
#endif
#ifdef RT_SOURCE_INJECTION
void rt_source_injection(void);
#endif

#ifdef RADTRANSFER
void rt_set_simple_inits(void);

#ifdef RT_DIFFUSION_CG
void rt_diffusion_cg_solve(void);
#endif

#ifdef RT_CHEM_PHOTOION
double rt_return_photon_number_density(MyIDType i, int k);
void rt_update_chemistry(void);
void rt_get_sigma(void);
double rt_GetCoolingTime(int i, double u, double rho);
double rt_cooling_photoheating(int i, double dt);
double rt_DoCooling(int, double);
double rt_DoHeating(int, double);
double rt_get_cooling_rate(int i, double entropy);
void rt_write_chemistry_stats(void);
void rt_get_lum_for_spectral_bin_stars(double T_eff, double luminosity_fraction[N_RT_FREQ_BINS]);
void rt_get_lum_gas(int target, double *je);
#endif

#endif


void find_block(char *label,FILE *fd);

#ifdef DISTORTIONTENSORPS
void get_half_kick_distortion(int pindex, MyBigFloat half_kick_add[6][6]);
void analyse_phase_space(int pindex,
                         MyBigFloat *s_1, MyBigFloat *s_2, MyBigFloat *s_3,
                         MyBigFloat *smear_x, MyBigFloat *smear_y, MyBigFloat *smear_z,
                         MyBigFloat *second_deriv, MyBigFloat *sigma);
MyBigFloat get_analytic_annihilation(MyBigFloat s_1, MyBigFloat s_2, MyBigFloat s_3,
                                   MyBigFloat s_1_prime, MyBigFloat s_2_prime, MyBigFloat s_3_prime,
                                   MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
                                   MyBigFloat dt, MyBigFloat sigma);
MyBigFloat get_max_caustic_density(MyBigFloat s_1, MyBigFloat s_2,
                                 MyBigFloat s_1_prime, MyBigFloat s_2_prime,
                                 MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
                                 MyBigFloat sigma,
                                 int pindex);
void get_current_ps_info(int pindex, MyBigFloat * flde, MyBigFloat * psde);
void do_phase_space_drift(int i, double dt_drift);
void do_the_phase_space_kick(int i, double dt_gravkick);
void do_long_range_phase_space_kick(int i, double dt_gravkick);
/* some math functions we need from phasespace_math.c */
void ludcmp(MyBigFloat ** a, int n, int *indx, MyBigFloat * d);
MyBigFloat **matrix(long nrl, long nrh, long ncl, long nch);
MyBigFloat *vector(long nl, long nh);
void free_matrix(MyBigFloat ** m, long nrl, long nrh, long ncl, long nch);
void free_vector(MyBigFloat * v, long nl, long nh);
void mult_matrix(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension, MyBigFloat ** matrix_result);
void mult_matrix_transpose_A(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension, MyBigFloat ** matrix_result);
void luinvert(MyBigFloat ** input_matrix, int n, MyBigFloat ** inverse_matrix);
void eigsrt(MyBigFloat d[], MyBigFloat ** v, int n);
void jacobi(MyBigFloat ** a, int n, MyBigFloat d[], MyBigFloat ** v, int *nrot);
#ifdef PMGRID
void pmtidaltensor_periodic_diff(void);
void pmtidaltensor_periodic_fourier(int component);
int pmtidaltensor_nonperiodic_diff(int grnr);
int pmtidaltensor_nonperiodic_fourier(int component, int grnr);
void check_tidaltensor_periodic(int particle_ID);
void check_tidaltensor_nonperiodic(int particle_ID);
#endif
#endif

#ifdef SCFPOTENTIAL
void SCF_do_center_of_mass_correction(double fac_rad, double start_rad, double fac_part, int max_iter);
void SCF_write(int task);
void SCF_calc_from_random(long *seed);
void SCF_calc_from_particles(void);
void SCF_init(void);
void SCF_reset(void);
void SCF_free(void);
void SCF_evaluate(MyDouble x, MyDouble y, MyDouble z, MyDouble *potential, MyDouble *ax, MyDouble *ay, MyDouble *az);
void SCF_collect_update(void);

void sphere_acc(double x, double y, double z, double *xa, double *ya, double *za);
void to_unit(double x, double y, double z, double *xs, double *ys, double *zs);
double ran1(long *idum);
double gasdev(long *idum);
double factrl(int n);
int nlm_all(int num, int n, int l, int m);
int nlm(int n, int l, int m);
int nl(int n, int l);
int lm(int l, int m);
int kdelta(int a, int b);
double gnlm_var(int n, int l, int m);
double hnlm_var(int n, int l, int m);
#endif

int ags_gravity_kernel_shared_check(short int particle_type_primary, short int particle_type_secondary);
#ifdef ADAPTIVE_GRAVSOFT_FORALL
void ags_setup_smoothinglengths(void);
int ags_ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
                                      int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                      int *ngblist, int type_of_searching_particle);
void ags_density(void);
int ags_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *ags_density_evaluate_primary(void *p);
void *ags_density_evaluate_secondary(void *p);
int ags_density_isactive(MyIDType i);
#endif

#ifdef ALTERNATIVE_PSORT
void init_sort_ID(MyIDType *data, int ndata);
#endif

void hydro_gradient_calc(void);
int GasGrad_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int gradient_iteration);
void *GasGrad_evaluate_primary(void *p, int gradient_iteration);
void *GasGrad_evaluate_secondary(void *p, int gradient_iteration);

#ifdef PARTICLE_EXCISION
void apply_excision();
#endif

#ifdef SIDM
#include "./sidm/sidm_proto.h"
#endif
