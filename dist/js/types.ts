/** Schema dist/js/schema/3pse/db/nist_jarvis/2024.3.13/atoms.json */

export interface NISTJARVISDbEntryAtomsKeySchemaBasedOnHttpsFigshareComArticlesDatasetMonolayerDataForHeterostructure22344571 {
  /**
   * Crystal lattice vectors as a 3x3 matrix, in Angstroms
   *
   * @minItems 3
   * @maxItems 3
   */
  lattice_mat?: [[number, number, number], [number, number, number], [number, number, number]];
  /**
   * Atomic coordinates for each atom in the unit cell
   *
   * @minItems 1
   */
  coords?: [[number, number, number], ...[number, number, number][]];
  /**
   * Atomic elements for each atom in the unit cell in the same order as `coords`
   *
   * @minItems 1
   */
  elements?: [string, ...string[]];
  /**
   * @minItems 3
   * @maxItems 3
   */
  abc?: [number, number, number];
  /**
   * @minItems 3
   * @maxItems 3
   */
  angles?: [number, number, number];
  /**
   * True if the coordinates are in Cartesian space, false if in fractional space
   */
  cartesian?: boolean;
  /**
   * Additional properties for each of the atoms
   */
  props?: string[];
}
 
/** Schema dist/js/schema/3pse/db/nist_jarvis/2024.3.13/db_entry.json */

export interface NISTJARVISDbEntrySchemaBasedOnHttpsFigshareComArticlesDatasetMonolayerDataForHeterostructure22344571 {
  atoms?: {
    /**
     * Crystal lattice vectors as a 3x3 matrix, in Angstroms
     *
     * @minItems 3
     * @maxItems 3
     */
    lattice_mat?: [[number, number, number], [number, number, number], [number, number, number]];
    /**
     * Atomic coordinates for each atom in the unit cell
     *
     * @minItems 1
     */
    coords?: [[number, number, number], ...[number, number, number][]];
    /**
     * Atomic elements for each atom in the unit cell in the same order as `coords`
     *
     * @minItems 1
     */
    elements?: [string, ...string[]];
    /**
     * @minItems 3
     * @maxItems 3
     */
    abc?: [number, number, number];
    /**
     * @minItems 3
     * @maxItems 3
     */
    angles?: [number, number, number];
    /**
     * True if the coordinates are in Cartesian space, false if in fractional space
     */
    cartesian?: boolean;
    /**
     * Additional properties for each of the atoms
     */
    props?: string[];
  };
  /**
   * The id of the entry in the database, e.g. JVASP-677
   */
  jid?: string;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/atomic_positions.json */

/**
 * https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1493
 */
export interface AtomicPositionsSchema {
  card_option?: "alat" | "bohr" | "angstrom" | "crystal" | "crystal_sg";
  values?: {
    /**
     * label of the atom as specified in ATOMIC_SPECIES
     */
    X?: string;
    /**
     * atomic positions
     */
    x: number;
    /**
     * atomic positions
     */
    y: number;
    /**
     * atomic positions
     */
    z: number;
    "if_pos(1)"?: number;
    "if_pos(2)"?: number;
    "if_pos(3)"?: number;
  }[];
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/atomic_species.json */

export interface AtomicSpeciesSchema {
  values?: {
    /**
     * label of the atom. Acceptable syntax: chemical symbol X (1 or 2 characters, case-insensitive) or chemical symbol plus a number or a letter, as in "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h; max total length cannot exceed 3 characters)
     */
    X?: string;
    /**
     * mass of the atomic species [amu: mass of C = 12]. Used only when performing Molecular Dynamics run or structural optimization runs using Damped MD. Not actually used in all other cases (but stored in data files, so phonon calculations will use these values unless other values are provided)
     */
    Mass_X?: number;
    /**
     * PseudoPot_X
     */
    PseudoPot_X?: string;
  }[];
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/cell.json */

export type CellSchema = CellSchema1 & CellSchema2;
export type CellSchema2 =
  | {
      /**
       * CASE ( calculation == 'vc-relax' )
       */
      cell_dynamics?: "none" | "sd" | "damp-pr" | "damp-w" | "bfgs";
    }
  | {
      /**
       * CASE ( calculation == 'vc-md' )
       */
      cell_dynamics?: "none" | "pr" | "w";
    };

export interface CellSchema1 {
  /**
   * Target pressure [KBar] in a variable-cell md or relaxation run.
   */
  press?: number;
  /**
   * Fictitious cell mass [amu] for variable-cell simulations (both 'vc-md' and 'vc-relax'). Default: 0.75*Tot_Mass/pi**2 for Parrinello-Rahman MD; 0.75*Tot_Mass/pi**2/Omega**(2/3) for Wentzcovitch MD
   */
  wmass?: number;
  /**
   * Used in the construction of the pseudopotential tables. It should exceed the maximum linear contraction of the cell during a simulation. Default: 2.0 for variable-cell calculations, 1.0 otherwise
   */
  cell_factor?: number;
  /**
   * Convergence threshold on the pressure for variable cell relaxation ('vc-relax' : note that the other convergence thresholds for ionic relaxation apply as well).
   */
  press_conv_thr?: number;
  /**
   * Select which of the cell parameters should be moved
   */
  cell_dofree?:
    | "all"
    | "ibrav"
    | "a"
    | "b"
    | "c"
    | "fixa"
    | "fixb"
    | "fixc"
    | "x"
    | "y"
    | "xy"
    | "xz"
    | "xyz"
    | "shape"
    | "volume"
    | "2Dxy"
    | "2Dshape"
    | "epitaxial_ab"
    | "epitaxial_ac"
    | "epitaxial_bc";
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/cell_parameters.json */

export interface CellParametersSchema {
  /**
   * label of the atom. Acceptable syntax: chemical symbol X (1 or 2 characters, case-insensitive) or chemical symbol plus a number or a letter, as in "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h; max total length cannot exceed 3 characters)
   */
  card_option?: "alat" | "bohr" | "angstrom";
  values?: {
    /**
     * @minItems 3
     * @maxItems 3
     */
    v1?: [number, number, number];
    /**
     * @minItems 3
     * @maxItems 3
     */
    v2?: [number, number, number];
    /**
     * @minItems 3
     * @maxItems 3
     */
    v3?: [number, number, number];
  };
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/control.json */

export interface ControlSchema {
  /**
   * A string describing the task to be performed
   */
  calculation?: "scf" | "nscf" | "bands" | "relax" | "md" | "vc-relax" | "vc-md";
  /**
   * Currently two verbosity levels are implemented: high, low. 'debug' and 'medium' have the same effect as 'high'; 'default' and 'minimal' as 'low'
   */
  verbosity?: "high" | "low" | "debug" | "medium" | "minimal" | "default";
  restart_mode?: "from_scratch" | "restart";
  /**
   * OBSOLETE - NO LONGER IMPLEMENTED
   */
  wf_collect?: boolean;
  /**
   * Default: 1 if calculation == 'scf', 'nscf', 'bands'; 50 for the other cases; Number of molecular-dynamics or structural optimization steps performed in this run. If set to 0, the code performs a quick "dry run", stopping just after initialization. This is useful to check for input correctness and to have the summary printed. NOTE: in MD calculations, the code will perform "nstep" steps even if restarting from a previously interrupted calculation.
   */
  nstep?: number;
  /**
   * band energies are written every iprint iterations
   */
  iprint?: number;
  /**
   * calculate stress. It is set to .TRUE. automatically if calculation == 'vc-md' or 'vc-relax'
   */
  tstress?: boolean;
  /**
   * calculate forces. It is set to .TRUE. automatically if calculation == 'relax','md','vc-md'
   */
  tprnfor?: boolean;
  /**
   * time step for molecular dynamics, in Rydberg atomic units (1 a.u.=4.8378 * 10^-17 s : beware, the CP code uses Hartree atomic units, half that much!!!)
   */
  dt?: number;
  /**
   * input, temporary, output files are found in this directory, see also wfcdir
   */
  outdir?: string;
  /**
   * This directory specifies where to store files generated by each processor (*.wfc{N}, *.igk{N}, etc.). Useful for machines without a parallel file system: set wfcdir to a local file system, while outdir should be a parallel or network file system, visible to all processors. Beware: in order to restart from interrupted runs, or to perform further calculations using the produced data files, you may need to copy files to outdir. Works only for pw.x.
   */
  wfcdir?: string;
  /**
   * prepended to input/output filenames: prefix.wfc, prefix.rho, etc.
   */
  prefix?: string;
  /**
   * OBSOLETE - NO LONGER IMPLEMENTED
   */
  lkpoint_dir?: boolean;
  /**
   * Jobs stops after max_seconds CPU time. Use this option in conjunction with option restart_mode if you need to split a job too long to complete into shorter jobs that fit into your batch queues.
   */
  max_seconds?: number;
  /**
   * Convergence threshold on total energy (a.u) for ionic minimization: the convergence criterion is satisfied when the total energy changes less than etot_conv_thr between two consecutive scf steps. Note that etot_conv_thr is extensive, like the total energy. See also forc_conv_thr - both criteria must be satisfied
   */
  etot_conv_thr?: number;
  /**
   * Convergence threshold on forces (a.u) for ionic minimization: the convergence criterion is satisfied when all components of all forces are smaller than forc_conv_thr. See also etot_conv_thr - both criteria must be satisfied
   */
  forc_conv_thr?: number;
  /**
   * Specifies the amount of disk I/O activity: (only for binary files and xml data file in data directory; other files printed at each molecular dynamics / structural optimization step are not controlled by this option )
   */
  disk_io?: "high" | "medium" | "low" | "nowf" | "none";
  /**
   * directory containing pseudopotential files. Default: value of the $ESPRESSO_PSEUDO environment variable if set; '$HOME/espresso/pseudo/' otherwise
   */
  pseudo_dir?: string;
  /**
   * If .TRUE. a saw-like potential simulating an electric field is added to the bare ionic potential. See variables edir, eamp, emaxpos, eopreg for the form and size of the added potential.
   */
  tefield?: boolean;
  /**
   * If .TRUE. and tefield==.TRUE. a dipole correction is also added to the bare ionic potential - implements the recipe of L. Bengtsson, PRB 59, 12301 (1999). See variables edir, emaxpos, eopreg for the form of the correction. Must be used ONLY in a slab geometry, for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE.
   */
  dipfield?: boolean;
  /**
   * If .TRUE. a homogeneous finite electric field described through the modern theory of the polarization is applied. This is different from tefield == .true. !
   */
  lelfield?: boolean;
  /**
   * In the case of a finite electric field  ( lelfield == .TRUE. ) it defines the number of iterations for converging the wavefunctions in the electric field Hamiltonian, for each external iteration on the charge density
   */
  nberrycyc?: number;
  /**
   * If .TRUE. perform orbital magnetization calculation.
   */
  lorbm?: boolean;
  /**
   * If .TRUE. perform a Berry phase calculation. See the header of PW/src/bp_c_phase.f90 for documentation
   */
  lberry?: boolean;
  /**
   * For Berry phase calculation: direction of the k-point strings in reciprocal space. Allowed values: 1, 2, 3 1=first, 2=second, 3=third reciprocal lattice vector For calculations with finite electric fields (lelfield==.true.) "gdir" is the direction of the field.
   */
  gdir?: number;
  /**
   * For Berry phase calculation: number of k-points to be calculated along each symmetry-reduced string. The same for calculation with finite electric fields (lelfield==.true.).
   */
  nppstr?: number;
  /**
   * In the case of charged cells (tot_charge .ne. 0) setting gate = .TRUE. represents the counter charge (i.e. -tot_charge) not by a homogeneous background charge but with a charged plate, which is placed at zgate (see below). Details of the gate potential can be found in T. Brumme, M. Calandra, F. Mauri; PRB 89, 245406 (2014). Note, that in systems which are not symmetric with respect to the plate, one needs to enable the dipole correction! (dipfield=.true.). Currently, symmetry can be used with gate=.true. but carefully check that no symmetry is included which maps z to -z even if in principle one could still use them for symmetric systems (i.e. no dipole correction). For nosym=.false. verbosity is set to 'high'. Note: this option was called "monopole" in v6.0 and 6.1 of pw.x
   */
  gate?: boolean;
  /**
   * IF .TRUE. , a two chemical potential calculation for the simulation of photoexcited systems is performed, constraining a fraction of the electrons in the conduction manifold.
   */
  twochem?: boolean;
  /**
   * If .TRUE. perform a constant bias potential (constant-mu) calculation for a system with ESM method. See the header of PW/src/fcp_module.f90 for documentation. To perform the calculation, you must set a namelist FCP.
   */
  lfcp?: boolean;
  /**
   * If .TRUE. perform a 3D-RISM-SCF calculation [for details see H.Sato et al., JCP 112, 9463 (2000), doi:10.1063/1.481564]. The solvent's distributions are calculated by 3D-RISM, though solute is treated as SCF. The charge density and the atomic positions are optimized, simultaneously with the solvents. To perform the calculation, you must set a namelist RISM and a card SOLVENTS. If assume_isolated = 'esm' and esm_bc = 'bc1', Laue-RISM is calculated instead of 3D-RISM and coupled with ESM method (i.e. ESM-RISM). [for details see S.Nishihara and M.Otani, PRB 96, 115429 (2017)]. The default of mixing_beta is 0.2 for both 3D-RISM and Laue-RISM. For structural relaxation with BFGS, ignore_wolfe is always .TRUE. .
   */
  trism?: boolean;
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/electrons.json */

export interface ElectronsSchema {
  /**
   * maximum number of iterations in a scf step. If exact exchange is active, this will affect the inner loops.
   */
  electron_maxstep?: number;
  /**
   * maximum number of outer iterations in a scf calculation with exact exchange.
   */
  exx_maxstep?: number;
  /**
   * If .false. do not stop molecular dynamics or ionic relaxation when electron_maxstep is reached. Use with care.
   */
  scf_must_converge?: boolean;
  conv_thr?: number;
  /**
   * If .TRUE. this turns on the use of an adaptive conv_thr for the inner scf loops when using EXX.
   */
  adaptive_thr?: boolean;
  /**
   * When adaptive_thr = .TRUE. this is the convergence threshold used for the first scf cycle.
   */
  conv_thr_init?: number;
  /**
   * When adaptive_thr = .TRUE. the convergence threshold for each scf cycle is given by: max( conv_thr, conv_thr_multi * dexx )
   */
  conv_thr_multi?: number;
  mixing_mode?: "plain" | "TF" | "local-TF";
  /**
   * mixing factor for self-consistency
   */
  mixing_beta?: number;
  /**
   * number of iterations used in mixing scheme
   */
  mixing_ndim?: number;
  /**
   * For DFT+U : number of iterations with fixed ns ( ns is the atomic density appearing in the Hubbard term ).
   */
  mixing_fixed_ns?: number;
  diagonalization?: "david" | "cg" | "ppcg" | "paro" | "ParO" | "rmm-davidson" | "rmm-paro";
  /**
   * Convergence threshold (ethr) for iterative diagonalization (the check is on eigenvalue convergence).
   */
  diago_thr_init?: number;
  /**
   * For conjugate gradient diagonalization:  max number of iterations
   */
  diago_cg_maxiter?: number;
  /**
   * For ppcg diagonalization:  max number of iterations
   */
  diago_ppcg_maxiter?: number;
  /**
   * For Davidson diagonalization: dimension of workspace (number of wavefunction packets, at least 2 needed).
   */
  diago_david_ndim?: number;
  /**
   * For RMM-DIIS diagonalization: dimension of workspace (number of wavefunction packets, at least 2 needed).
   */
  diago_rmm_ndim?: number;
  /**
   * If .TRUE., RMM-DIIS is performed up to converge. If .FALSE., RMM-DIIS is performed only once.
   */
  diago_rmm_conv?: boolean;
  /**
   * For RMM-DIIS diagonalization: blocking size of Gram-Schmidt orthogonalization
   */
  diago_gs_nblock?: number;
  /**
   * If .TRUE. all the empty states are diagonalized at the same level of accuracy of the occupied ones. Otherwise the empty states are diagonalized using a larger threshold (this should not affect total energy, forces, and other ground-state properties).
   */
  diago_full_acc?: boolean;
  /**
   * Amplitude of the finite electric field (in Ry a.u.; 1 a.u. = 36.3609*10^10 V/m). Used only if lelfield==.TRUE. and if k-points (K_POINTS card) are not automatic.
   */
  efield?: number;
  /**
   * @minItems 3
   * @maxItems 3
   */
  efield_cart?: [number, number, number];
  efield_phase?: "read" | "write" | "none";
  startingpot?: "atomic" | "file";
  startingwfc?: "atomic" | "atomic+random" | "random" | "file";
  /**
   * If .true., use a real-space algorithm for augmentation charges of ultrasoft pseudopotentials and PAWsets. Faster but numerically less accurate than the default G-space algorithm. Use with care and after testing!
   */
  tqr?: boolean;
  /**
   * If .true., exploit real-space localization to compute matrix elements for nonlocal projectors. Faster and in principle better scaling than the default G-space algorithm, but numerically less accurate, may lead to some loss of translational invariance. Use with care and after testing!
   */
  real_space?: boolean;
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/hubbard.json */

export interface HubbardSchema {
  card_option?: "atomic" | "ortho-atomic" | "norm-atomic" | "wf" | "pseudo";
  values?:
    | (
        | {
            /**
             * string constant "U"; indicates the specs for the U parameter will be given
             */
            U?: "U";
            /**
             * label of the atom (as defined in ATOMIC_SPECIES)
             */
            label?: string;
            /**
             * specs of the manifold (e.g., 3d, 2p...)
             */
            manifold?: string;
            /**
             * value of the U parameter (in eV)
             */
            u_val?: number;
          }
        | {
            /**
             * string constant "J0"; indicates the specs for the J0 parameter will be given
             */
            J0?: "J0";
            /**
             * label of the atom (as defined in ATOMIC_SPECIES)
             */
            label?: string;
            /**
             * specs of the manifold (e.g., 3d, 2p...)
             */
            manifold?: string;
            /**
             * value of the J0 parameter (in eV)
             */
            j0_val?: number;
          }
      )[]
    | {
        /**
         * character describing the type of Hubbard parameter allowed values: U, J and either B (for d-orbitals) or E2 and E3 (for f-orbitals)
         */
        paramType?: "U" | "J" | "B" | "E2" | "E3";
        /**
         * label of the atom (as defined in ATOMIC_SPECIES)
         */
        label?: string;
        /**
         * specs of the manifold (e.g., 3d, 2p...)
         */
        manifold?: string;
        /**
         * value of the J0 parameter (in eV)
         */
        paramValue?: number;
      }[]
    | (
        | {
            /**
             * string constant "U"; indicates the specs for the U parameter will be given
             */
            U?: "U";
            /**
             * label of the atom (as defined in ATOMIC_SPECIES)
             */
            label?: string;
            /**
             * specs of the manifold (e.g., 3d, 2p...)
             */
            manifold?: string;
            /**
             * value of the U parameter (in eV)
             */
            u_val?: number;
          }
        | {
            /**
             * string constant "J0"; indicates the specs for the J0 parameter will be given
             */
            J0?: "J0";
            /**
             * label of the atom (as defined in ATOMIC_SPECIES)
             */
            label?: string;
            /**
             * specs of the manifold (e.g., 3d, 2p...)
             */
            manifold?: string;
            /**
             * value of the J0 parameter (in eV)
             */
            j0_val?: number;
          }
        | {
            /**
             * string constant "V"; indicates the specs for the V parameter will be given
             */
            V?: "V";
            /**
             * label of the atom I (as defined in ATOMIC_SPECIES)
             */
            "label(I)"?: string;
            /**
             * specs of the manifold for atom I (e.g., 3d, 2p...)
             */
            "manifold(I)"?: string;
            /**
             * label of the atom J (as defined in ATOMIC_SPECIES)
             */
            "label(J)"?: string;
            /**
             * specs of the manifold for atom J (e.g., 3d, 2p...)
             */
            "manifold(J)"?: string;
            /**
             * index of the atom I
             */
            I?: number;
            /**
             * index of the atom J
             */
            J?: number;
            /**
             * value of the V parameter for the atom pair I,J (in eV)
             */
            "v_val(I,J)"?: number;
          }
      )[];
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/ions.json */

export type IonsSchema = IonsSchema1 & IonsSchema2;
export type IonsSchema2 =
  | {
      /**
       * CASE: calculation == 'relax'
       */
      ion_dynamics?: "bfgs" | "damp" | "fire";
    }
  | {
      /**
       * CASE: calculation == 'md'
       */
      ion_dynamics?: "verlet" | "langevin" | "langevin-smc";
    }
  | {
      /**
       * CASE: calculation == 'vc-relax'
       */
      ion_dynamics?: "bfgs" | "damp";
    }
  | {
      /**
       * CASE: calculation == 'vc-md'
       */
      ion_dynamics?: "beeman";
    };

export interface IonsSchema1 {
  ion_positions?: "default" | "from_input";
  ion_velocities?: "default" | "from_input";
  /**
   * Used to extrapolate the potential from preceding ionic steps.
   */
  pot_extrapolation?: "none" | "atomic" | "first_order" | "second_order";
  /**
   * Used to extrapolate the wavefunctions from preceding ionic steps.
   */
  wfc_extrapolation?: "none" | "first_order" | "second_order";
  /**
   * This keyword is useful when simulating the dynamics and/or the thermodynamics of an isolated system. If set to true the total torque of the internal forces is set to zero by adding new forces that compensate the spurious interaction with the periodic images. This allows for the use of smaller supercells.
   */
  remove_rigid_rot?: boolean;
  ion_temperature?:
    | "rescaling"
    | "rescale-v"
    | "rescale-T"
    | "reduce-T"
    | "berendsen"
    | "andersen"
    | "svr"
    | "initial"
    | "not_controlled";
  /**
   * Starting temperature (Kelvin) in MD runs target temperature for most thermostats.
   */
  tempw?: number;
  /**
   * Tolerance for velocity rescaling. Velocities are rescaled if the run-averaged and target temperature differ more than tolp.
   */
  tolp?: number;
  delta_t?: number;
  nraise?: number;
  /**
   * This keyword applies only in the case of molecular dynamics or damped dynamics. If true the ions are refolded at each step into the supercell.
   */
  refold_pos?: boolean;
  /**
   * Max reduction factor for conv_thr during structural optimization conv_thr is automatically reduced when the relaxation approaches convergence so that forces are still accurate, but conv_thr will not be reduced to less that conv_thr / upscale.
   */
  upscale?: number;
  /**
   * Number of old forces and displacements vectors used in the PULAY mixing of the residual vectors obtained on the basis of the inverse hessian matrix given by the BFGS algorithm.
   */
  bfgs_ndim?: number;
  /**
   * Maximum ionic displacement in the structural relaxation. (bfgs only)
   */
  trust_radius_max?: number;
  /**
   * Minimum ionic displacement in the structural relaxation BFGS is reset when trust_radius < trust_radius_min. (bfgs only)
   */
  trust_radius_min?: number;
  /**
   * Initial ionic displacement in the structural relaxation. (bfgs only)
   */
  trust_radius_ini?: number;
  w_1?: number;
  /**
   * Parameters used in line search based on the Wolfe conditions. (bfgs only)
   */
  w_2?: number;
  /**
   * Initial value of the alpha mixing factor in the FIRE minimization scheme; recommended values are between 0.1 and 0.3
   */
  fire_alpha_init?: number;
  /**
   * Scaling of the alpha mixing parameter for steps with P > 0;
   */
  fire_falpha?: number;
  /**
   * Minimum number of steps with P > 0 before increase of dt
   */
  fire_nmin?: number;
  /**
   * Factor for increasing dt
   */
  fire_f_inc?: number;
  /**
   * Factor for decreasing dt
   */
  fire_f_dec?: number;
  /**
   * Determines the maximum value of dt in the FIRE minimization; dtmax = fire_dtmax*dt
   */
  fire_dtmax?: number;
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/k_points.json */

export interface KPointsSchema {
  card_option?: "tpiba" | "automatic" | "crystal" | "gamma" | "tpiba_b" | "crystal_b" | "tpiba_c" | "crystal_c";
  values?:
    | {
        /**
         * Number of supplied special k-points.
         */
        nks?: number;
        xk_x?: number;
        xk_y?: number;
        xk_z?: number;
        wk?: number;
      }[]
    | {
        /**
         * Number of supplied special k-points.
         */
        nk1?: number;
        /**
         * Number of supplied special k-points.
         */
        nk2?: number;
        /**
         * Number of supplied special k-points.
         */
        nk3?: number;
        /**
         * Number of supplied special k-points.
         */
        sk1?: number;
        /**
         * Number of supplied special k-points.
         */
        sk2?: number;
        /**
         * Number of supplied special k-points.
         */
        sk3?: number;
      }
    | null;
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x/system.json */

export type SystemSchema = SystemSchema1 & SystemSchema2;
export type SystemSchema1 =
  | {
      /**
       * @minItems 6
       * @maxItems 6
       */
      celldm?: [number, number, number, number, number, number];
    }
  | {
      A?: number;
      B?: number;
      C?: number;
      cosAB?: number;
      cosAC?: number;
      cosBC?: number;
    };

export interface SystemSchema2 {
  ibrav: number;
  /**
   * number of atoms in the unit cell (ALL atoms, except if space_group is set, in which case, INEQUIVALENT atoms)
   */
  nat: number;
  /**
   * number of types of atoms in the unit cell
   */
  ntyp: number;
  /**
   * Default: for an insulator, nbnd = number of valence bands (nbnd = # of electrons /2); for a metal, 20% more (minimum 4 more)
   */
  nbnd?: number;
  /**
   * Default: nbnd_cond = nbnd - # of electrons / 2 in the collinear case; nbnd_cond = nbnd - # of electrons in the noncollinear case.
   */
  nbnd_cond?: number;
  tot_charge?: number;
  /**
   * starting charge on atomic type 'i', to create starting potential with startingpot = 'atomic'.
   */
  starting_charge?: number;
  /**
   * Total majority spin charge - minority spin charge. Used to impose a specific total electronic magnetization. If unspecified then tot_magnetization variable is ignored and the amount of electronic magnetization is determined during the self-consistent cycle.
   */
  tot_magnetization?: number;
  starting_magnetization?: number[];
  /**
   * kinetic energy cutoff (Ry) for wavefunctions
   */
  ecutwfc: number;
  /**
   * Kinetic energy cutoff (Ry) for charge density and potential For norm-conserving pseudopotential you should stick to the default value, you can reduce it by a little but it will introduce noise especially on forces and stress. Default: 4 * ecutwfc
   */
  ecutrho?: number;
  /**
   * Kinetic energy cutoff (Ry) for the exact exchange operator in EXX type calculations. By default this is the same as ecutrho but in some EXX calculations, a significant speed-up can be obtained by reducing ecutfock, at the expense of some loss in accuracy. Must be .gt. ecutwfc. Not implemented for stress calculation and for US-PP and PAW pseudopotentials.
   */
  ecutfock?: number;
  /**
   * Three-dimensional FFT mesh (hard grid) for charge density (and scf potential). If not specified the grid is calculated based on the cutoff for charge density (see also ecutrho)
   */
  nr1?: number;
  /**
   * Three-dimensional FFT mesh (hard grid) for charge density (and scf potential). If not specified the grid is calculated based on the cutoff for charge density (see also ecutrho)
   */
  nr2?: number;
  /**
   * Three-dimensional FFT mesh (hard grid) for charge density (and scf potential). If not specified the grid is calculated based on the cutoff for charge density (see also ecutrho)
   */
  nr3?: number;
  /**
   * Three-dimensional mesh for wavefunction FFT and for the smooth part of charge density ( smooth grid ). Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
   */
  nr1s?: number;
  /**
   * Three-dimensional mesh for wavefunction FFT and for the smooth part of charge density ( smooth grid ). Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
   */
  nr2s?: number;
  /**
   * Three-dimensional mesh for wavefunction FFT and for the smooth part of charge density ( smooth grid ). Coincides with nr1, nr2, nr3 if ecutrho = 4 * ecutwfc ( default )
   */
  nr3s?: number;
  nosym?: boolean;
  nosym_evc?: boolean;
  /**
   * if (.TRUE.) disable the usage of k => -k symmetry (time reversal) in k-point generation
   */
  noinv?: boolean;
  /**
   * if (.TRUE.) disable the usage of magnetic symmetry operations that consist in a rotation + time reversal.
   */
  no_t_rev?: boolean;
  /**
   * if (.TRUE.) force the symmetry group to be symmorphic by disabling symmetry operations having an associated fractionary translation
   */
  force_symmorphic?: boolean;
  use_all_frac?: boolean;
  occupations?: "smearing" | "tetrahedra" | "tetrahedra_lin" | "tetrahedra_opt" | "fixed" | "from_input";
  one_atom_occupations?: boolean;
  starting_spin_angle?: boolean;
  /**
   * value of the gaussian spreading (Ry) for brillouin-zone integration in the conduction manifold in a two-chemical potential calculation (twochem=.true.).
   */
  degauss_cond?: number;
  /**
   * Number of electrons placed in the conduction manifold in a two-chemical potential calculation (twochem=.true.). Of the total # of electrons nelec, nelec-nelec_cond will occupy the valence manifold and nelec_cond will be constrained in the conduction manifold.
   */
  nelec_cond?: number;
  /**
   * value of the gaussian spreading (Ry) for brillouin-zone integration in metals.
   */
  degauss?: number;
  smearing?:
    | "gaussian"
    | "gauss"
    | "methfessel-paxton"
    | "m-p"
    | "mp"
    | "marzari-vanderbilt"
    | "cold"
    | "m-v"
    | "mv"
    | "fermi-dirac"
    | "f-d"
    | "fd";
  nspin?: number;
  /**
   * Strength of the gammaDFT potential.
   */
  sic_gamma?: number;
  /**
   * Type of polaron in gammaDFT.
   */
  pol_type?: "e" | "h";
  /**
   * Enable the calculation of the total energy in gammaDFT. When .true., a preliminary calculation is performed to calculate the electron density in the absence of the polaron. When .false., the total energy printed in output should not be considered. For structural relaxations, it is recommended to use .false. to avoid doubling the computational cost.
   */
  sic_energy?: boolean;
  /**
   * Valence band shift (in eV) through self-consistent scissor operator. When performing gammaDFT calculations of polarons, the polaron level is not shifted.
   */
  sci_vb?: number;
  /**
   * Conduction band band shift (in eV) through self-consistent scissor operator. When performing gammaDFT calculations of polarons, the polaron level is not shifted.
   */
  sci_cb?: number;
  /**
   * if .true. the program will perform a noncollinear calculation.
   */
  noncolin?: boolean;
  ecfixed?: number;
  qcutz?: number;
  q2sigma?: number;
  /**
   * Exchange-correlation functional: eg 'PBE', 'BLYP' etc See Modules/funct.f90 for allowed values. Overrides the value read from pseudopotential files. Use with care and if you know what you are doing!
   */
  input_dft?: string;
  /**
   * Use Adaptively Compressed Exchange operator as in Lin Lin, J. Chem. Theory Comput. 2016, 12, 2242--2249, doi:10.1021/acs.jctc.6b00092
   */
  ace?: boolean;
  /**
   * Fraction of EXX for hybrid functional calculations. In the case of input_dft='PBE0', the default value is 0.25, while for input_dft='B3LYP' the exx_fraction default value is 0.20.
   */
  exx_fraction?: number;
  /**
   * screening_parameter for HSE like hybrid functionals.
   */
  screening_parameter?: number;
  exxdiv_treatment?: "gygi-baldereschi" | "vcut_spherical" | "vcut_ws" | "none";
  /**
   * Specific for EXX. If .true., extrapolate the G=0 term of the potential
   */
  x_gamma_extrapolation?: boolean;
  /**
   * Reciprocal space cutoff for correcting Coulomb potential divergencies at small q vectors.
   */
  ecutvcut?: number;
  /**
   * Three-dimensional mesh for q (k1-k2) sampling of the Fock operator (EXX). Can be smaller than the number of k-points.
   */
  nqx1?: number;
  /**
   * Three-dimensional mesh for q (k1-k2) sampling of the Fock operator (EXX). Can be smaller than the number of k-points.
   */
  nqx2?: number;
  /**
   * Three-dimensional mesh for q (k1-k2) sampling of the Fock operator (EXX). Can be smaller than the number of k-points.
   */
  nqx3?: number;
  /**
   * Overlap threshold over which the exchange integral over a pair of localized orbitals is included in the evaluation of EXX operator. Any value greater than 0.0 triggers the SCDM localization and the evaluation on EXX using the localized orbitals. Very small value of the threshold should yield the same result as the default EXX evaluation
   */
  localization_thr?: number;
  Hubbard_occ?: [number, number, number][];
  Hubbard_alpha?: number[];
  Hubbard_beta?: number[];
  starting_ns_eigenvalue?: number[][][];
  /**
   * If true, nscf calculation will exit in restart mode, scf calculation will restart from there if DMFT updates are provided as hdf5 archive. Scf calculation should be used only with electron_maxstep = 1.
   */
  dmft?: boolean;
  /**
   * prepended to hdf5 archive: dmft_prefix.h5
   */
  dmft_prefix?: string;
  /**
   * If ensemble_energies = .true., an ensemble of xc energies is calculated non-selfconsistently for perturbed exchange-enhancement factors and LDA vs. PBE correlation ratios after each converged electronic ground state calculation.
   */
  ensemble_energies?: boolean;
  /**
   * The direction of the electric field or dipole correction is parallel to the bg(:,edir) reciprocal lattice vector, so the potential is constant in planes defined by FFT grid points; edir = 1, 2 or 3. Used only if tefield is .TRUE.
   */
  edir?: number;
  /**
   * Position of the maximum of the saw-like potential along crystal axis edir, within the  unit cell (see below), 0 < emaxpos < 1 Used only if tefield is .TRUE.
   */
  emaxpos?: number;
  /**
   * Zone in the unit cell where the saw-like potential decreases. ( see below, 0 < eopreg < 1 ). Used only if tefield is .TRUE.
   */
  eopreg?: number;
  eamp?: number;
  /**
   * The angle expressed in degrees between the initial magnetization and the z-axis. For noncollinear calculations only; index i runs over the atom types.
   *
   * @minItems 1
   * @maxItems 1
   */
  angle1?: [number];
  /**
   * The angle expressed in degrees between the projection of the initial magnetization on x-y plane and the x-axis. For noncollinear calculations only.
   *
   * @minItems 1
   * @maxItems 1
   */
  angle2?: [number];
  /**
   * When starting a non collinear calculation using an existing density file from a collinear lsda calculation assumes previous density points in z direction and rotates it in the direction described by angle1 and angle2 variables for atomic type 1
   */
  lforcet?: boolean;
  /**
   * Used to perform constrained calculations in magnetic systems.
   */
  constrained_magnetization?: "none" | "total" | "atomic" | "total direction" | "atomic direction";
  /**
   * @minItems 3
   * @maxItems 3
   */
  fixed_magnetization?: [number, number, number];
  /**
   * parameter used for constrained_magnetization calculations N.B.: if the scf calculation does not converge, try to reduce lambda to obtain convergence, then restart the run with a larger lambda
   */
  lambda?: number;
  /**
   * determines when atomic magnetic moments are printed on output
   */
  report?: number;
  /**
   * if .TRUE. the noncollinear code can use a pseudopotential with spin-orbit.
   */
  lspinorb?: boolean;
  /**
   * Used to perform calculation assuming the system to be isolated (a molecule or a cluster in a 3D supercell)
   */
  assume_isolated?: "none" | "makov-payne" | "m-p" | "mp" | "martyna-tuckerman" | "m-t" | "mt" | "esm" | "2D";
  /**
   * If assume_isolated = 'esm', determines the boundary conditions used for either side of the slab.
   */
  esm_bc?: "pbc" | "bc1" | "bc2" | "bc3";
  /**
   * If assume_isolated = 'esm', determines the position offset [in a.u.] of the start of the effective screening region, measured relative to the cell edge. (ESM region begins at z = +/- [L_z/2 + esm_w] ).
   */
  esm_w?: number;
  /**
   * If assume_isolated = 'esm' and esm_bc = 'bc2', gives the magnitude of the electric field [Ry/a.u.] to be applied between semi-infinite ESM electrodes.
   */
  esm_efield?: number;
  /**
   * If assume_isolated = 'esm', gives the number of z-grid points for the polynomial fit along the cell edge.
   */
  esm_nfit?: number;
  /**
   * If .TRUE. perform a constant bias potential (constant-mu) calculation with Grand-Canonical SCF.
   */
  lgcscf?: boolean;
  /**
   * The target Fermi energy (eV) of GC-SCF. One can start with appropriate total charge of the system by giving tot_charge
   */
  gcscf_mu?: number;
  /**
   * Convergence threshold of Fermi energy (eV) for GC-SCF.
   */
  gcscf_conv_thr?: number;
  /**
   * Mixing factor for GC-SCF. Larger values are recommended, if systems with small DOS on Fermi surface as graphite.
   */
  gcscf_beta?: number;
  /**
   * Type of Van der Waals correction
   */
  vdw_corr?:
    | "none"
    | "grimme-d2"
    | "Grimme-D2"
    | "DFT-D"
    | "dft-d"
    | "grimme-d3"
    | "Grimme-D3"
    | "DFT-D3"
    | "dft-d3"
    | "TS"
    | "ts"
    | "ts-vdw"
    | "ts-vdW"
    | "tkatchenko-scheffler"
    | "MBD"
    | "mbd"
    | "many-body-dispersion"
    | "mbd_vdw"
    | "XDM"
    | "xdm";
  /**
   * OBSOLESCENT, same as vdw_corr='DFT-D'
   */
  london?: boolean;
  /**
   * global scaling parameter for DFT-D. Default is good for PBE.
   */
  london_s6?: number;
  /**
   * atomic C6 coefficient of each atom type
   */
  london_c6?: number;
  /**
   * atomic vdw radii of each atom type
   */
  london_rvdw?: number;
  /**
   * cutoff radius (a.u.) for dispersion interactions
   */
  london_rcut?: number;
  /**
   * Version of Grimme implementation of Grimme-D3
   */
  dftd3_version?: number;
  /**
   * Turn three-body terms in Grimme-D3 on. If .false. two-body contributions only are computed, using two-body parameters of Grimme-D3. If dftd3_version=2, three-body contribution is always disabled.
   */
  dftd3_threebody?: boolean;
  /**
   * Optional: controls the convergence of the vdW energy (and forces). The default value is a safe choice, likely too safe, but you do not gain much in increasing it
   */
  ts_vdw_econv_thr?: number;
  /**
   * Optional: set it to .TRUE. when computing the Tkatchenko-Scheffler vdW energy or the Many-Body dispersion (MBD) energy for an isolated (non-periodic) system.
   */
  ts_vdw_isolated?: boolean;
  /**
   * OBSOLESCENT, same as vdw_corr='xdm'
   */
  xdm?: boolean;
  /**
   * Damping function parameter a1 (adimensional)
   */
  xdm_a1?: number;
  /**
   * Damping function parameter a2 (angstrom)
   */
  xdm_a2?: number;
  /**
   * The number of the space group of the crystal, as given in the International Tables of Crystallography A (ITA)
   */
  space_group?: number;
  /**
   * Used only for monoclinic lattices
   */
  uniqueb?: boolean;
  /**
   * Used only for space groups that in the ITA allow the use of two different origins
   */
  origin_choice?: number;
  /**
   * Used only for rhombohedral space groups.
   */
  rhombohedral?: boolean;
  /**
   * used only if gate = .TRUE.
   */
  zgate?: number;
  /**
   * used only if gate = .TRUE.
   */
  relaxz?: boolean;
  /**
   * used only if gate = .TRUE.
   */
  block?: boolean;
  /**
   * used only if gate = .TRUE. and block = .TRUE.
   */
  block_1?: number;
  /**
   * used only if gate = .TRUE. and block = .TRUE.
   */
  block_2?: number;
  /**
   * used only if gate = .TRUE. and block = .TRUE.
   */
  block_height?: number;
  /**
   * Number of activated external ionic force fields.
   */
  nextffield?: number;
}
 
/** Schema dist/js/schema/3pse/file/applications/espresso/7.2/pw.x.json */

export interface PwxMainSchema {
  "&CONTROL"?: {
    /**
     * A string describing the task to be performed
     */
    calculation?: "scf" | "nscf" | "bands" | "relax" | "md" | "vc-relax" | "vc-md";
    /**
     * Currently two verbosity levels are implemented: high, low. 'debug' and 'medium' have the same effect as 'high'; 'default' and 'minimal' as 'low'
     */
    verbosity?: "high" | "low" | "debug" | "medium" | "minimal" | "default";
    restart_mode?: "from_scratch" | "restart";
    /**
     * OBSOLETE - NO LONGER IMPLEMENTED
     */
    wf_collect?: boolean;
    /**
     * Default: 1 if calculation == 'scf', 'nscf', 'bands'; 50 for the other cases; Number of molecular-dynamics or structural optimization steps performed in this run. If set to 0, the code performs a quick "dry run", stopping just after initialization. This is useful to check for input correctness and to have the summary printed. NOTE: in MD calculations, the code will perform "nstep" steps even if restarting from a previously interrupted calculation.
     */
    nstep?: number;
    /**
     * band energies are written every iprint iterations
     */
    iprint?: number;
    /**
     * calculate stress. It is set to .TRUE. automatically if calculation == 'vc-md' or 'vc-relax'
     */
    tstress?: boolean;
    /**
     * calculate forces. It is set to .TRUE. automatically if calculation == 'relax','md','vc-md'
     */
    tprnfor?: boolean;
    /**
     * time step for molecular dynamics, in Rydberg atomic units (1 a.u.=4.8378 * 10^-17 s : beware, the CP code uses Hartree atomic units, half that much!!!)
     */
    dt?: number;
    /**
     * input, temporary, output files are found in this directory, see also wfcdir
     */
    outdir?: string;
    /**
     * This directory specifies where to store files generated by each processor (*.wfc{N}, *.igk{N}, etc.). Useful for machines without a parallel file system: set wfcdir to a local file system, while outdir should be a parallel or network file system, visible to all processors. Beware: in order to restart from interrupted runs, or to perform further calculations using the produced data files, you may need to copy files to outdir. Works only for pw.x.
     */
    wfcdir?: string;
    /**
     * prepended to input/output filenames: prefix.wfc, prefix.rho, etc.
     */
    prefix?: string;
    /**
     * OBSOLETE - NO LONGER IMPLEMENTED
     */
    lkpoint_dir?: boolean;
    /**
     * Jobs stops after max_seconds CPU time. Use this option in conjunction with option restart_mode if you need to split a job too long to complete into shorter jobs that fit into your batch queues.
     */
    max_seconds?: number;
    /**
     * Convergence threshold on total energy (a.u) for ionic minimization: the convergence criterion is satisfied when the total energy changes less than etot_conv_thr between two consecutive scf steps. Note that etot_conv_thr is extensive, like the total energy. See also forc_conv_thr - both criteria must be satisfied
     */
    etot_conv_thr?: number;
    /**
     * Convergence threshold on forces (a.u) for ionic minimization: the convergence criterion is satisfied when all components of all forces are smaller than forc_conv_thr. See also etot_conv_thr - both criteria must be satisfied
     */
    forc_conv_thr?: number;
    /**
     * Specifies the amount of disk I/O activity: (only for binary files and xml data file in data directory; other files printed at each molecular dynamics / structural optimization step are not controlled by this option )
     */
    disk_io?: "high" | "medium" | "low" | "nowf" | "none";
    /**
     * directory containing pseudopotential files. Default: value of the $ESPRESSO_PSEUDO environment variable if set; '$HOME/espresso/pseudo/' otherwise
     */
    pseudo_dir?: string;
    /**
     * If .TRUE. a saw-like potential simulating an electric field is added to the bare ionic potential. See variables edir, eamp, emaxpos, eopreg for the form and size of the added potential.
     */
    tefield?: boolean;
    /**
     * If .TRUE. and tefield==.TRUE. a dipole correction is also added to the bare ionic potential - implements the recipe of L. Bengtsson, PRB 59, 12301 (1999). See variables edir, emaxpos, eopreg for the form of the correction. Must be used ONLY in a slab geometry, for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE.
     */
    dipfield?: boolean;
    /**
     * If .TRUE. a homogeneous finite electric field described through the modern theory of the polarization is applied. This is different from tefield == .true. !
     */
    lelfield?: boolean;
    /**
     * In the case of a finite electric field  ( lelfield == .TRUE. ) it defines the number of iterations for converging the wavefunctions in the electric field Hamiltonian, for each external iteration on the charge density
     */
    nberrycyc?: number;
    /**
     * If .TRUE. perform orbital magnetization calculation.
     */
    lorbm?: boolean;
    /**
     * If .TRUE. perform a Berry phase calculation. See the header of PW/src/bp_c_phase.f90 for documentation
     */
    lberry?: boolean;
    /**
     * For Berry phase calculation: direction of the k-point strings in reciprocal space. Allowed values: 1, 2, 3 1=first, 2=second, 3=third reciprocal lattice vector For calculations with finite electric fields (lelfield==.true.) "gdir" is the direction of the field.
     */
    gdir?: number;
    /**
     * For Berry phase calculation: number of k-points to be calculated along each symmetry-reduced string. The same for calculation with finite electric fields (lelfield==.true.).
     */
    nppstr?: number;
    /**
     * In the case of charged cells (tot_charge .ne. 0) setting gate = .TRUE. represents the counter charge (i.e. -tot_charge) not by a homogeneous background charge but with a charged plate, which is placed at zgate (see below). Details of the gate potential can be found in T. Brumme, M. Calandra, F. Mauri; PRB 89, 245406 (2014). Note, that in systems which are not symmetric with respect to the plate, one needs to enable the dipole correction! (dipfield=.true.). Currently, symmetry can be used with gate=.true. but carefully check that no symmetry is included which maps z to -z even if in principle one could still use them for symmetric systems (i.e. no dipole correction). For nosym=.false. verbosity is set to 'high'. Note: this option was called "monopole" in v6.0 and 6.1 of pw.x
     */
    gate?: boolean;
    /**
     * IF .TRUE. , a two chemical potential calculation for the simulation of photoexcited systems is performed, constraining a fraction of the electrons in the conduction manifold.
     */
    twochem?: boolean;
    /**
     * If .TRUE. perform a constant bias potential (constant-mu) calculation for a system with ESM method. See the header of PW/src/fcp_module.f90 for documentation. To perform the calculation, you must set a namelist FCP.
     */
    lfcp?: boolean;
    /**
     * If .TRUE. perform a 3D-RISM-SCF calculation [for details see H.Sato et al., JCP 112, 9463 (2000), doi:10.1063/1.481564]. The solvent's distributions are calculated by 3D-RISM, though solute is treated as SCF. The charge density and the atomic positions are optimized, simultaneously with the solvents. To perform the calculation, you must set a namelist RISM and a card SOLVENTS. If assume_isolated = 'esm' and esm_bc = 'bc1', Laue-RISM is calculated instead of 3D-RISM and coupled with ESM method (i.e. ESM-RISM). [for details see S.Nishihara and M.Otani, PRB 96, 115429 (2017)]. The default of mixing_beta is 0.2 for both 3D-RISM and Laue-RISM. For structural relaxation with BFGS, ignore_wolfe is always .TRUE. .
     */
    trism?: boolean;
  };
  "&SYSTEM"?:
    | {
        /**
         * @minItems 6
         * @maxItems 6
         */
        celldm?: [number, number, number, number, number, number];
      }
    | {
        A?: number;
        B?: number;
        C?: number;
        cosAB?: number;
        cosAC?: number;
        cosBC?: number;
      };
  "&ELECTRONS"?: {
    /**
     * maximum number of iterations in a scf step. If exact exchange is active, this will affect the inner loops.
     */
    electron_maxstep?: number;
    /**
     * maximum number of outer iterations in a scf calculation with exact exchange.
     */
    exx_maxstep?: number;
    /**
     * If .false. do not stop molecular dynamics or ionic relaxation when electron_maxstep is reached. Use with care.
     */
    scf_must_converge?: boolean;
    conv_thr?: number;
    /**
     * If .TRUE. this turns on the use of an adaptive conv_thr for the inner scf loops when using EXX.
     */
    adaptive_thr?: boolean;
    /**
     * When adaptive_thr = .TRUE. this is the convergence threshold used for the first scf cycle.
     */
    conv_thr_init?: number;
    /**
     * When adaptive_thr = .TRUE. the convergence threshold for each scf cycle is given by: max( conv_thr, conv_thr_multi * dexx )
     */
    conv_thr_multi?: number;
    mixing_mode?: "plain" | "TF" | "local-TF";
    /**
     * mixing factor for self-consistency
     */
    mixing_beta?: number;
    /**
     * number of iterations used in mixing scheme
     */
    mixing_ndim?: number;
    /**
     * For DFT+U : number of iterations with fixed ns ( ns is the atomic density appearing in the Hubbard term ).
     */
    mixing_fixed_ns?: number;
    diagonalization?: "david" | "cg" | "ppcg" | "paro" | "ParO" | "rmm-davidson" | "rmm-paro";
    /**
     * Convergence threshold (ethr) for iterative diagonalization (the check is on eigenvalue convergence).
     */
    diago_thr_init?: number;
    /**
     * For conjugate gradient diagonalization:  max number of iterations
     */
    diago_cg_maxiter?: number;
    /**
     * For ppcg diagonalization:  max number of iterations
     */
    diago_ppcg_maxiter?: number;
    /**
     * For Davidson diagonalization: dimension of workspace (number of wavefunction packets, at least 2 needed).
     */
    diago_david_ndim?: number;
    /**
     * For RMM-DIIS diagonalization: dimension of workspace (number of wavefunction packets, at least 2 needed).
     */
    diago_rmm_ndim?: number;
    /**
     * If .TRUE., RMM-DIIS is performed up to converge. If .FALSE., RMM-DIIS is performed only once.
     */
    diago_rmm_conv?: boolean;
    /**
     * For RMM-DIIS diagonalization: blocking size of Gram-Schmidt orthogonalization
     */
    diago_gs_nblock?: number;
    /**
     * If .TRUE. all the empty states are diagonalized at the same level of accuracy of the occupied ones. Otherwise the empty states are diagonalized using a larger threshold (this should not affect total energy, forces, and other ground-state properties).
     */
    diago_full_acc?: boolean;
    /**
     * Amplitude of the finite electric field (in Ry a.u.; 1 a.u. = 36.3609*10^10 V/m). Used only if lelfield==.TRUE. and if k-points (K_POINTS card) are not automatic.
     */
    efield?: number;
    /**
     * @minItems 3
     * @maxItems 3
     */
    efield_cart?: [number, number, number];
    efield_phase?: "read" | "write" | "none";
    startingpot?: "atomic" | "file";
    startingwfc?: "atomic" | "atomic+random" | "random" | "file";
    /**
     * If .true., use a real-space algorithm for augmentation charges of ultrasoft pseudopotentials and PAWsets. Faster but numerically less accurate than the default G-space algorithm. Use with care and after testing!
     */
    tqr?: boolean;
    /**
     * If .true., exploit real-space localization to compute matrix elements for nonlocal projectors. Faster and in principle better scaling than the default G-space algorithm, but numerically less accurate, may lead to some loss of translational invariance. Use with care and after testing!
     */
    real_space?: boolean;
  };
  "&IONS"?:
    | (
        | {
            /**
             * CASE: calculation == 'relax'
             */
            ion_dynamics?: "bfgs" | "damp" | "fire";
          }
        | {
            /**
             * CASE: calculation == 'md'
             */
            ion_dynamics?: "verlet" | "langevin" | "langevin-smc";
          }
        | {
            /**
             * CASE: calculation == 'vc-relax'
             */
            ion_dynamics?: "bfgs" | "damp";
          }
        | {
            /**
             * CASE: calculation == 'vc-md'
             */
            ion_dynamics?: "beeman";
          }
      )
    | null;
  "&CELL"?:
    | (
        | {
            /**
             * CASE ( calculation == 'vc-relax' )
             */
            cell_dynamics?: "none" | "sd" | "damp-pr" | "damp-w" | "bfgs";
          }
        | {
            /**
             * CASE ( calculation == 'vc-md' )
             */
            cell_dynamics?: "none" | "pr" | "w";
          }
      )
    | null;
  ATOMIC_SPECIES?: {
    values?: {
      /**
       * label of the atom. Acceptable syntax: chemical symbol X (1 or 2 characters, case-insensitive) or chemical symbol plus a number or a letter, as in "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h; max total length cannot exceed 3 characters)
       */
      X?: string;
      /**
       * mass of the atomic species [amu: mass of C = 12]. Used only when performing Molecular Dynamics run or structural optimization runs using Damped MD. Not actually used in all other cases (but stored in data files, so phonon calculations will use these values unless other values are provided)
       */
      Mass_X?: number;
      /**
       * PseudoPot_X
       */
      PseudoPot_X?: string;
    }[];
  };
  /**
   * https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1493
   */
  ATOMIC_POSITIONS?: {
    card_option?: "alat" | "bohr" | "angstrom" | "crystal" | "crystal_sg";
    values?: {
      /**
       * label of the atom as specified in ATOMIC_SPECIES
       */
      X?: string;
      /**
       * atomic positions
       */
      x: number;
      /**
       * atomic positions
       */
      y: number;
      /**
       * atomic positions
       */
      z: number;
      "if_pos(1)"?: number;
      "if_pos(2)"?: number;
      "if_pos(3)"?: number;
    }[];
  };
  K_POINTS?: {
    card_option?: "tpiba" | "automatic" | "crystal" | "gamma" | "tpiba_b" | "crystal_b" | "tpiba_c" | "crystal_c";
    values?:
      | {
          /**
           * Number of supplied special k-points.
           */
          nks?: number;
          xk_x?: number;
          xk_y?: number;
          xk_z?: number;
          wk?: number;
        }[]
      | {
          /**
           * Number of supplied special k-points.
           */
          nk1?: number;
          /**
           * Number of supplied special k-points.
           */
          nk2?: number;
          /**
           * Number of supplied special k-points.
           */
          nk3?: number;
          /**
           * Number of supplied special k-points.
           */
          sk1?: number;
          /**
           * Number of supplied special k-points.
           */
          sk2?: number;
          /**
           * Number of supplied special k-points.
           */
          sk3?: number;
        }
      | null;
  };
  CELL_PARAMETERS?: {
    /**
     * label of the atom. Acceptable syntax: chemical symbol X (1 or 2 characters, case-insensitive) or chemical symbol plus a number or a letter, as in "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h; max total length cannot exceed 3 characters)
     */
    card_option?: "alat" | "bohr" | "angstrom";
    values?: {
      /**
       * @minItems 3
       * @maxItems 3
       */
      v1?: [number, number, number];
      /**
       * @minItems 3
       * @maxItems 3
       */
      v2?: [number, number, number];
      /**
       * @minItems 3
       * @maxItems 3
       */
      v3?: [number, number, number];
    };
  };
  HUBBARD?: {
    card_option?: "atomic" | "ortho-atomic" | "norm-atomic" | "wf" | "pseudo";
    values?:
      | (
          | {
              /**
               * string constant "U"; indicates the specs for the U parameter will be given
               */
              U?: "U";
              /**
               * label of the atom (as defined in ATOMIC_SPECIES)
               */
              label?: string;
              /**
               * specs of the manifold (e.g., 3d, 2p...)
               */
              manifold?: string;
              /**
               * value of the U parameter (in eV)
               */
              u_val?: number;
            }
          | {
              /**
               * string constant "J0"; indicates the specs for the J0 parameter will be given
               */
              J0?: "J0";
              /**
               * label of the atom (as defined in ATOMIC_SPECIES)
               */
              label?: string;
              /**
               * specs of the manifold (e.g., 3d, 2p...)
               */
              manifold?: string;
              /**
               * value of the J0 parameter (in eV)
               */
              j0_val?: number;
            }
        )[]
      | {
          /**
           * character describing the type of Hubbard parameter allowed values: U, J and either B (for d-orbitals) or E2 and E3 (for f-orbitals)
           */
          paramType?: "U" | "J" | "B" | "E2" | "E3";
          /**
           * label of the atom (as defined in ATOMIC_SPECIES)
           */
          label?: string;
          /**
           * specs of the manifold (e.g., 3d, 2p...)
           */
          manifold?: string;
          /**
           * value of the J0 parameter (in eV)
           */
          paramValue?: number;
        }[]
      | (
          | {
              /**
               * string constant "U"; indicates the specs for the U parameter will be given
               */
              U?: "U";
              /**
               * label of the atom (as defined in ATOMIC_SPECIES)
               */
              label?: string;
              /**
               * specs of the manifold (e.g., 3d, 2p...)
               */
              manifold?: string;
              /**
               * value of the U parameter (in eV)
               */
              u_val?: number;
            }
          | {
              /**
               * string constant "J0"; indicates the specs for the J0 parameter will be given
               */
              J0?: "J0";
              /**
               * label of the atom (as defined in ATOMIC_SPECIES)
               */
              label?: string;
              /**
               * specs of the manifold (e.g., 3d, 2p...)
               */
              manifold?: string;
              /**
               * value of the J0 parameter (in eV)
               */
              j0_val?: number;
            }
          | {
              /**
               * string constant "V"; indicates the specs for the V parameter will be given
               */
              V?: "V";
              /**
               * label of the atom I (as defined in ATOMIC_SPECIES)
               */
              "label(I)"?: string;
              /**
               * specs of the manifold for atom I (e.g., 3d, 2p...)
               */
              "manifold(I)"?: string;
              /**
               * label of the atom J (as defined in ATOMIC_SPECIES)
               */
              "label(J)"?: string;
              /**
               * specs of the manifold for atom J (e.g., 3d, 2p...)
               */
              "manifold(J)"?: string;
              /**
               * index of the atom I
               */
              I?: number;
              /**
               * index of the atom J
               */
              J?: number;
              /**
               * value of the V parameter for the atom pair I,J (in eV)
               */
              "v_val(I,J)"?: number;
            }
        )[];
  };
}
 
/** Schema dist/js/schema/core/abstract/2d_data.json */

export interface DimensionDataSchema {
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/core/abstract/2d_plot.json */

export interface DimensionPlotSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: string;
    /**
     * units for an axis
     */
    units?: string;
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: string;
    /**
     * units for an axis
     */
    units?: string;
  };
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/core/abstract/3d_grid.json */

export interface DimensionalGridSchema {
  /**
   * @minItems 3
   * @maxItems 3
   */
  dimensions: [number, number, number];
  /**
   * @minItems 3
   * @maxItems 3
   */
  shifts: [number, number, number];
}
 
/** Schema dist/js/schema/core/abstract/3d_tensor.json */

/**
 * @minItems 3
 * @maxItems 3
 */
export type DimensionalTensorSchema = [[number, number, number], [number, number, number], [number, number, number]];
 
/** Schema dist/js/schema/core/abstract/3d_vector_basis.json */

export interface DimensionalVectorBasis {
  /**
   * @minItems 3
   * @maxItems 3
   */
  a: [number, number, number];
  /**
   * @minItems 3
   * @maxItems 3
   */
  b: [number, number, number];
  /**
   * @minItems 3
   * @maxItems 3
   */
  c: [number, number, number];
}
 
/** Schema dist/js/schema/core/abstract/point.json */

/**
 * @minItems 3
 * @maxItems 3
 */
export type PointSchema = [number, number, number];
 
/** Schema dist/js/schema/core/abstract/vector.json */

export type VectorSchema = [number, number, number] | [boolean, boolean, boolean];
 
/** Schema dist/js/schema/core/primitive/1d_data_series.json */

export type DimensionDataSeriesSchema = [number | string, ...(number | string)[]][];
 
/** Schema dist/js/schema/core/primitive/3d_lattice.json */

export interface DimensionalLatticeSchema {
  /**
   * length of the first lattice vector
   */
  a: number;
  /**
   * length of the second lattice vector
   */
  b: number;
  /**
   * length of the third lattice vector
   */
  c: number;
  /**
   * angle between first and second lattice vector
   */
  alpha: number;
  /**
   * angle between second and third lattice vector
   */
  beta: number;
  /**
   * angle between first and third lattice vector
   */
  gamma: number;
}
 
/** Schema dist/js/schema/core/primitive/array_of_3_booleans.json */

/**
 * @minItems 3
 * @maxItems 3
 */
export type ArrayOf3BooleanElementsSchema = [boolean, boolean, boolean];
 
/** Schema dist/js/schema/core/primitive/array_of_3_numbers.json */

/**
 * @minItems 3
 * @maxItems 3
 */
export type ArrayOf3NumberElementsSchema = [number, number, number];
 
/** Schema dist/js/schema/core/primitive/array_of_ids.json */

/**
 * array of objects containing integer id each
 */
export type AtomicIds = {
  /**
   * integer id of this entry
   */
  id?: number;
}[];
 
/** Schema dist/js/schema/core/primitive/array_of_strings.json */

/**
 * array of strings, e.g. metadata tags
 */
export type ArrayOfStrings = string[];
 
/** Schema dist/js/schema/core/primitive/axis.json */

export interface AxisSchema {
  /**
   * label of an axis object
   */
  label: string;
  /**
   * units for an axis
   */
  units?: string;
}
 
/** Schema dist/js/schema/core/primitive/group_info.json */

export interface GroupInfoSchemaForNodesInAGraph {
  /**
   * Human-readable name of group of nodes
   */
  groupName?: string;
  /**
   * Unique identifier of the group a node belongs to
   */
  groupId?: string;
}
 
/** Schema dist/js/schema/core/primitive/integer_one_or_zero.json */

export type IntegerOneOrZero = number;
 
/** Schema dist/js/schema/core/primitive/linked_list/base_node.json */

export interface BasicNodeSchemaLinkedList {
  /**
   * Flowchart ID of next node
   */
  next?: string;
  /**
   * Whether node is head node or not
   */
  head?: boolean;
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
}
/**
 * This interface was referenced by `BasicNodeSchemaLinkedList`'s JSON-Schema
 * via the `definition` "flowchart-id".
 */
export interface FlowchartId {
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
}
 
/** Schema dist/js/schema/core/primitive/linked_list/named_node.json */

export interface NamedNodeSchema {
  /**
   * Flowchart ID of next node
   */
  next?: string;
  /**
   * Whether node is head node or not
   */
  head?: boolean;
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
  /**
   * entity name
   */
  name?: string;
}
/**
 * This interface was referenced by `NamedNodeSchema`'s JSON-Schema
 * via the `definition` "flowchart-id".
 */
export interface FlowchartId {
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
}
 
/** Schema dist/js/schema/core/primitive/linked_list/named_node_in_group.json */

export interface NamedNodeInGroupSchema {
  /**
   * Flowchart ID of next node
   */
  next?: string;
  /**
   * Whether node is head node or not
   */
  head?: boolean;
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Human-readable name of group of nodes
   */
  groupName?: string;
  /**
   * Unique identifier of the group a node belongs to
   */
  groupId?: string;
}
/**
 * This interface was referenced by `NamedNodeInGroupSchema`'s JSON-Schema
 * via the `definition` "flowchart-id".
 */
export interface FlowchartId {
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
}
 
/** Schema dist/js/schema/core/primitive/linked_list/node_with_type.json */

export interface TypedNodeSchema {
  type?: string;
  /**
   * Flowchart ID of next node
   */
  next?: string;
  /**
   * Whether node is head node or not
   */
  head?: boolean;
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
}
/**
 * This interface was referenced by `TypedNodeSchema`'s JSON-Schema
 * via the `definition` "flowchart-id".
 */
export interface FlowchartId {
  /**
   * Unique flowchart ID of node
   */
  flowchartId: string;
}
 
/** Schema dist/js/schema/core/primitive/linked_list.json */

export type LinkedListSchema = (
  | {
      /**
       * Flowchart ID of next node
       */
      next?: string;
      /**
       * Whether node is head node or not
       */
      head?: boolean;
      /**
       * Unique flowchart ID of node
       */
      flowchartId: string;
    }
  | {
      /**
       * Flowchart ID of next node
       */
      next?: string;
      /**
       * Whether node is head node or not
       */
      head?: boolean;
      /**
       * Unique flowchart ID of node
       */
      flowchartId: string;
      /**
       * entity name
       */
      name?: string;
    }
  | {
      /**
       * Flowchart ID of next node
       */
      next?: string;
      /**
       * Whether node is head node or not
       */
      head?: boolean;
      /**
       * Unique flowchart ID of node
       */
      flowchartId: string;
      /**
       * entity name
       */
      name?: string;
      /**
       * Human-readable name of group of nodes
       */
      groupName?: string;
      /**
       * Unique identifier of the group a node belongs to
       */
      groupId?: string;
    }
  | {
      type?: string;
      /**
       * Flowchart ID of next node
       */
      next?: string;
      /**
       * Whether node is head node or not
       */
      head?: boolean;
      /**
       * Unique flowchart ID of node
       */
      flowchartId: string;
    }
)[];
 
/** Schema dist/js/schema/core/primitive/scalar.json */

export interface ScalarSchema {
  value: number;
}
 
/** Schema dist/js/schema/core/primitive/slugified_entry.json */

/**
 * container for machine- and human-readable identifier
 */
export interface SlugifiedEntry {
  /**
   * descriptive human-readable name of entry
   */
  name: string;
  /**
   * machine-readable identifier
   */
  slug: string;
}
 
/** Schema dist/js/schema/core/primitive/slugified_entry_or_slug.json */

/**
 * contains either object with slugified entry or slug only as a string
 */
export type SlugifiedEntryOrSlug =
  | {
      /**
       * descriptive human-readable name of entry
       */
      name: string;
      /**
       * machine-readable identifier
       */
      slug: string;
    }
  | string;
 
/** Schema dist/js/schema/core/primitive/string.json */

export interface PrimitiveString {
  value: string;
}
 
/** Schema dist/js/schema/core/reference/exabyte.json */

export interface CoreReferenceExabyte {
  /**
   * Material's identity. Used for protoProperties.
   */
  materialId?: string;
  /**
   * Job's identity
   */
  jobId?: string;
  /**
   * Id of the unit that extracted the result
   */
  unitId?: string;
}
 
/** Schema dist/js/schema/core/reference/experiment/condition.json */

export interface ConditionSchema {
  /**
   * condition unit
   */
  units?: string;
  /**
   * array of condition values
   */
  scalar?: {
    value?: string;
  }[];
  /**
   * human-readable name of the condition
   */
  name: string;
}
 
/** Schema dist/js/schema/core/reference/experiment/location.json */

export interface LocationSchema {
  /**
   * location latitude
   */
  latitude: number;
  /**
   * location longitude
   */
  longitude: number;
}
 
/** Schema dist/js/schema/core/reference/experiment.json */

export interface InfoForCharacteristicObtainedByExperiment {
  type?: "experiment";
  /**
   * experiment authors
   */
  authors: {
    first: string;
    middle?: string;
    last: string;
    affiliation?: string;
  }[];
  /**
   * method used in experiment
   */
  method: string;
  conditions: {
    /**
     * condition unit
     */
    units?: string;
    /**
     * array of condition values
     */
    scalar?: {
      value?: string;
    }[];
    /**
     * human-readable name of the condition
     */
    name: string;
  }[];
  location?: {
    /**
     * location latitude
     */
    latitude: number;
    /**
     * location longitude
     */
    longitude: number;
  };
  /**
   * epoch time.
   */
  timestamp: number;
  /**
   * Note about experiment
   */
  note?: string;
  /**
   * references to literature articles
   */
  references?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  }[];
}
 
/** Schema dist/js/schema/core/reference/literature/name.json */

export interface ExperimentAuthorSchema {
  first: string;
  middle?: string;
  last: string;
  affiliation?: string;
}
 
/** Schema dist/js/schema/core/reference/literature/pages.json */

export interface PagesSchema {
  start: string;
  end?: string;
}
 
/** Schema dist/js/schema/core/reference/literature.json */

export interface LiteratureReferenceSchema {
  type?: "literature";
  /**
   * Digital Object Identifier of the reference.
   */
  doi?: string;
  /**
   * International Standard Book Number of the reference.
   */
  isbn?: string;
  /**
   * International Standard Serial Number of the reference.
   */
  issn?: string;
  /**
   * Internet address of the reference.
   */
  url?: string;
  /**
   * Publisher of the work.
   */
  publisher?: string;
  /**
   * Journal in which the work appeared.
   */
  journal?: string;
  /**
   * Volume of the series in which the work appeared.
   */
  volume?: string;
  /**
   * Year in which the reference was published.
   */
  year?: string;
  /**
   * Issue of the collection in which the work appeared.
   */
  issue?: string;
  /**
   * Start and end pages of the work.
   */
  pages?: {
    start: string;
    end?: string;
  };
  /**
   * List of authors of the work.
   */
  authors?: {
    first: string;
    middle?: string;
    last: string;
    affiliation?: string;
  }[];
  /**
   * List of editors of the work.
   */
  editors?: {
    first: string;
    middle?: string;
    last: string;
    affiliation?: string;
  }[];
  /**
   * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
   */
  reference?: {}[];
}
 
/** Schema dist/js/schema/core/reference/modeling/exabyte.json */

export interface InfoForCharacteristicObtainedByExabyteCalculation {
  type?: "exabyte";
  /**
   * job identifier
   */
  _id: string;
  owner: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
}
 
/** Schema dist/js/schema/core/reference/modeling.json */

export type InfoForPropertyObtainedByModelingOnlySupportsExabyteOriginatedDataAtmButEasilyExtendable = {
  type?: "exabyte";
  /**
   * job identifier
   */
  _id: string;
  owner: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
};
 
/** Schema dist/js/schema/core/reference.json */

export type ReferenceSchemaUsingAnyOfInsteadOfOneOfBelowBCCurrentReferenceSchemasOverlap =
  | {
      type?: "exabyte";
      /**
       * job identifier
       */
      _id: string;
      owner: {
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity class
         */
        cls?: string;
        /**
         * entity slug
         */
        slug?: string;
      };
    }
  | {
      type?: "experiment";
      /**
       * experiment authors
       */
      authors: {
        first: string;
        middle?: string;
        last: string;
        affiliation?: string;
      }[];
      /**
       * method used in experiment
       */
      method: string;
      conditions: {
        /**
         * condition unit
         */
        units?: string;
        /**
         * array of condition values
         */
        scalar?: {
          value?: string;
        }[];
        /**
         * human-readable name of the condition
         */
        name: string;
      }[];
      location?: {
        /**
         * location latitude
         */
        latitude: number;
        /**
         * location longitude
         */
        longitude: number;
      };
      /**
       * epoch time.
       */
      timestamp: number;
      /**
       * Note about experiment
       */
      note?: string;
      /**
       * references to literature articles
       */
      references?: {
        type?: "literature";
        /**
         * Digital Object Identifier of the reference.
         */
        doi?: string;
        /**
         * International Standard Book Number of the reference.
         */
        isbn?: string;
        /**
         * International Standard Serial Number of the reference.
         */
        issn?: string;
        /**
         * Internet address of the reference.
         */
        url?: string;
        /**
         * Publisher of the work.
         */
        publisher?: string;
        /**
         * Journal in which the work appeared.
         */
        journal?: string;
        /**
         * Volume of the series in which the work appeared.
         */
        volume?: string;
        /**
         * Year in which the reference was published.
         */
        year?: string;
        /**
         * Issue of the collection in which the work appeared.
         */
        issue?: string;
        /**
         * Start and end pages of the work.
         */
        pages?: {
          start: string;
          end?: string;
        };
        /**
         * List of authors of the work.
         */
        authors?: {
          first: string;
          middle?: string;
          last: string;
          affiliation?: string;
        }[];
        /**
         * List of editors of the work.
         */
        editors?: {
          first: string;
          middle?: string;
          last: string;
          affiliation?: string;
        }[];
        /**
         * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
         */
        reference?: {}[];
      }[];
    }
  | {
      type?: "literature";
      /**
       * Digital Object Identifier of the reference.
       */
      doi?: string;
      /**
       * International Standard Book Number of the reference.
       */
      isbn?: string;
      /**
       * International Standard Serial Number of the reference.
       */
      issn?: string;
      /**
       * Internet address of the reference.
       */
      url?: string;
      /**
       * Publisher of the work.
       */
      publisher?: string;
      /**
       * Journal in which the work appeared.
       */
      journal?: string;
      /**
       * Volume of the series in which the work appeared.
       */
      volume?: string;
      /**
       * Year in which the reference was published.
       */
      year?: string;
      /**
       * Issue of the collection in which the work appeared.
       */
      issue?: string;
      /**
       * Start and end pages of the work.
       */
      pages?: {
        start: string;
        end?: string;
      };
      /**
       * List of authors of the work.
       */
      authors?: {
        first: string;
        middle?: string;
        last: string;
        affiliation?: string;
      }[];
      /**
       * List of editors of the work.
       */
      editors?: {
        first: string;
        middle?: string;
        last: string;
        affiliation?: string;
      }[];
      /**
       * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
       */
      reference?: {}[];
    };
 
/** Schema dist/js/schema/core/reusable/atomic_data/per_orbital.json */

/**
 * Atomic properties per orbital e.g., Hubbard U parameters.
 */
export interface AtomicDataPerOrbital {
  /**
   * Site number or index in the lattice
   */
  id?: number;
  /**
   * Example: Co1, Mn
   */
  atomicSpecies?: string;
  orbitalName?: string;
}
 
/** Schema dist/js/schema/core/reusable/atomic_data/per_orbital_pair.json */

/**
 * Atomic properties per orbital pair e.g., Hubbard V parameters.
 */
export interface AtomicDataPerOrbitalPair {
  /**
   * Site number or index in the lattice
   */
  id?: number;
  /**
   * Site number or index in the lattice of second site
   */
  id2?: number;
  /**
   * Example: Co1, Mn
   */
  atomicSpecies?: string;
  /**
   * Example: Co2, O
   */
  atomicSpecies2?: string;
  orbitalName?: string;
  orbitalName2?: string;
  /**
   * Distance between two sites in Bohr.
   */
  distance?: number;
}
 
/** Schema dist/js/schema/core/reusable/atomic_data/value_number.json */

/**
 * Numeric value specific to atomic data
 */
export interface AtomicDataNumericProperties {
  /**
   * Value related to a specific property, e.g., Hubbard U, V etc.
   */
  value?: number;
}
 
/** Schema dist/js/schema/core/reusable/atomic_data/value_string.json */

/**
 * String value specific to atomic data
 */
export interface AtomicDataStringProperties {
  /**
   * String value specific to atomic data
   */
  value?: string;
}
 
/** Schema dist/js/schema/core/reusable/atomic_data_per_orbital_numeric.json */

/**
 * Atomic properties per orbital pair with numeric value e.g., Hubbard V parameters.
 */
export interface AtomicDataPerOrbitalNumeric {
  /**
   * Site number or index in the lattice
   */
  id?: number;
  /**
   * Example: Co1, Mn
   */
  atomicSpecies?: string;
  orbitalName?: string;
  /**
   * Value related to a specific property, e.g., Hubbard U, V etc.
   */
  value?: number;
}
 
/** Schema dist/js/schema/core/reusable/atomic_data_per_orbital_pair_numeric.json */

/**
 * Atomic properties per orbital pair with numeric value e.g., Hubbard V parameters.
 */
export interface AtomicDataPerOrbitalPairNumeric {
  /**
   * Site number or index in the lattice
   */
  id?: number;
  /**
   * Site number or index in the lattice of second site
   */
  id2?: number;
  /**
   * Example: Co1, Mn
   */
  atomicSpecies?: string;
  /**
   * Example: Co2, O
   */
  atomicSpecies2?: string;
  orbitalName?: string;
  orbitalName2?: string;
  /**
   * Distance between two sites in Bohr.
   */
  distance?: number;
  /**
   * Value related to a specific property, e.g., Hubbard U, V etc.
   */
  value?: number;
}
 
/** Schema dist/js/schema/core/reusable/atomic_orbital.json */

export interface AtomicOrbitalSchema {
  orbitalName?: string;
  orbitalIndex?: number;
  principalNumber?: number;
  angularMomentum?: number;
  /**
   * Shell occupation
   */
  occupation?: number;
}
 
/** Schema dist/js/schema/core/reusable/atomic_scalars.json */

/**
 * array of objects containing integer id each
 */
export type AtomicScalarsVectorsSchema = {
  value?: {
    value: number;
  };
  /**
   * integer id of this entry
   */
  id?: number;
}[];
 
/** Schema dist/js/schema/core/reusable/atomic_strings.json */

/**
 * array of objects containing integer id each
 */
export type AtomicStringsVectorsSchema = {
  value?: string;
  /**
   * integer id of this entry
   */
  id?: number;
}[];
 
/** Schema dist/js/schema/core/reusable/atomic_vectors.json */

/**
 * array of objects containing integer id each
 */
export type AtomicVectorsSchema = {
  value?: [number, number, number] | [boolean, boolean, boolean];
  /**
   * integer id of this entry
   */
  id?: number;
}[];
 
/** Schema dist/js/schema/core/reusable/band_gap.json */

export interface BandGapSchema {
  /**
   * @minItems 3
   * @maxItems 3
   */
  kpointConduction?: [number, number, number];
  /**
   * @minItems 3
   * @maxItems 3
   */
  kpointValence?: [number, number, number];
  /**
   * eigenvalue at k-point in conduction band
   */
  eigenvalueConduction?: number;
  /**
   * eigenvalue at k-point in valence band
   */
  eigenvalueValence?: number;
  spin?: number;
  type: "direct" | "indirect";
  units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
  value: number;
}
 
/** Schema dist/js/schema/core/reusable/categories.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface ReusableCategoriesSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/core/reusable/category_path.json */

/**
 * TODO: Use regex once schema draft version has been updated
 */
export type CategoryPathSchema = string;
 
/** Schema dist/js/schema/core/reusable/dielectric_tensor_component.json */

/**
 * Schema for a function of frequency yielding a nx3 matrix
 */
export interface DielectricTensor {
  /**
   * Real or imaginary part of the dielectric tensor component
   */
  part: "real" | "imaginary";
  spin?: number;
  /**
   * Frequencies
   */
  frequencies: number[];
  /**
   * Matrix with 3 columns, e.g. x, y, z
   */
  components: [number, number, number][];
}
 
/** Schema dist/js/schema/core/reusable/energy.json */

export interface EnergySchema {
  name: string;
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/core/reusable/file_metadata.json */

export interface FileMetadata {
  /**
   * Relative path to the directory that contains the file.
   */
  pathname?: string;
  /**
   * Basename of the file
   */
  basename?: string;
  /**
   * What kind of file this is, e.g. image / text
   */
  filetype?: string;
}
 
/** Schema dist/js/schema/core/reusable/frequency_function_matrix.json */

/**
 * Schema for a function of frequency yielding a nx3 matrix
 */
export interface CoreReusableFrequencyFunctionMatrix {
  /**
   * Frequencies
   */
  frequencies?: number[];
  /**
   * Matrix with 3 columns, e.g. x, y, z
   */
  components?: [number, number, number][];
}
 
/** Schema dist/js/schema/core/reusable/object_storage_container_data.json */

export interface ObjectStorageContainerData {
  /**
   * Object storage container for the file
   */
  CONTAINER?: string;
  /**
   * Name of the file inside the object storage bucket
   */
  NAME?: string;
  /**
   * Object storage provider
   */
  PROVIDER?: string;
  /**
   * Region for the object container specified in Container
   */
  REGION?: string;
  /**
   * Size of the file in bytes
   */
  SIZE?: number;
  /**
   * Unix timestamp showing when the file was last modified
   */
  TIMESTAMP?: string;
}
 
/** Schema dist/js/schema/definitions/units.json */

export interface DefinitionsUnits {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/element.json */

export interface ElementSchema {
  /**
   * Element symbol.
   */
  symbol?: string;
  /**
   * list of elemental properties
   */
  properties?: (
    | {
        name?: "atomic_radius";
        units?:
          | "km"
          | "m"
          | "pm"
          | "nm"
          | "angstrom"
          | "a.u."
          | "bohr"
          | "fractional"
          | "crystal"
          | "cartesian"
          | "alat";
        value: number;
      }
    | {
        name?: "electronegativity";
        value: number;
      }
    | {
        name?: "ionization_potential";
        units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
        value: number;
      }
  )[];
}
 
/** Schema dist/js/schema/in_memory_entity/base.json */

export interface BaseInMemoryEntitySchema {
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
}
 
/** Schema dist/js/schema/in_memory_entity/defaultable.json */

export interface DefaultableInMemoryEntitySchema {
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
}
 
/** Schema dist/js/schema/in_memory_entity/named.json */

export interface NamedInMemoryEntitySchema {
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
}
 
/** Schema dist/js/schema/in_memory_entity/named_defaultable.json */

export interface NamedDefaultableInMemoryEntitySchema {
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
}
 
/** Schema dist/js/schema/in_memory_entity/named_defaultable_has_metadata.json */

export interface NamedDefaultableHasMetadataInMemoryEntitySchema {
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  metadata?: {};
}
 
/** Schema dist/js/schema/in_memory_entity/named_defaultable_runtime_items.json */

export interface NamedDefaultableRuntimeItemsInMemoryEntitySchema {
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
}
 
/** Schema dist/js/schema/job/base.json */

export interface JobBaseSchema {
  /**
   * Identity used to track jobs originated from command-line
   */
  rmsId?: string;
  /**
   * job status
   */
  status:
    | "pre-submission"
    | "queued"
    | "submitted"
    | "active"
    | "finished"
    | "terminate-queued"
    | "terminated"
    | "error"
    | "deleted"
    | "timeout";
  /**
   * Approximate start time of the job. e.g. within 10 min
   */
  startTime?: string;
  /**
   * The path to the working directory of this job, when the job originates from command-line
   */
  workDir?: string;
  /**
   * Custom keywords prefixed with validate correspond to custom validation methods implemented downstream
   */
  compute: {
    /**
     * Name of the submission queues: https://docs.mat3ra.com/infrastructure/resource/queues/. Below enums are for Azure, then AWS circa 2022-08, hence the duplication.
     */
    queue:
      | "D"
      | "OR"
      | "OF"
      | "OFplus"
      | "SR"
      | "SF"
      | "SFplus"
      | "GPOF"
      | "GP2OF"
      | "GP4OF"
      | "GPSF"
      | "GP2SF"
      | "GP4SF"
      | "OR4"
      | "OR8"
      | "OR16"
      | "SR4"
      | "SR8"
      | "SR16"
      | "GOF"
      | "G4OF"
      | "G8OF"
      | "GSF"
      | "G4SF"
      | "G8SF";
    /**
     * number of nodes used for the job inside the RMS.
     */
    nodes: number;
    /**
     * number of CPUs used for the job inside the RMS.
     */
    ppn: number;
    /**
     * Wallclock time limit for computing a job. Clock format: 'hh:mm:ss'
     */
    timeLimit: string;
    /**
     * Convention to use when reasoning about time limits
     */
    timeLimitType?: "per single attempt" | "compound";
    /**
     * Job is allowed to restart on termination.
     */
    isRestartable?: boolean;
    /**
     * Email notification for the job: n - never, a - job aborted, b - job begins, e - job ends. Last three could be combined.
     */
    notify?: string;
    /**
     * Email address to notify about job execution.
     */
    email?: string;
    /**
     * Maximum CPU count per node. This parameter is used to let backend job submission infrastructure know that this job is to be charged for the maximum CPU per node instead of the actual ppn. For premium/fast queues where resources are provisioned on-demand and exclusively per user.
     */
    maxCPU?: number;
    /**
     * Optional arguments specific to using application - VASP, Quantum Espresso, etc. Specified elsewhere
     */
    arguments?: {
      /**
       * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
       */
      nimage?: number;
      /**
       * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
       */
      npools?: number;
      /**
       * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
       */
      nband?: number;
      /**
       * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
       */
      ntg?: number;
      /**
       * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
       */
      ndiag?: number;
    };
    /**
     * Cluster where the job is executed. Optional on create. Required on job submission.
     */
    cluster?: {
      /**
       * FQDN of the cluster. e.g. master-1-staging.exabyte.io
       */
      fqdn?: string;
      /**
       * Job's identity in RMS. e.g. 1234.master-1-staging.exabyte.io
       */
      jid?: string;
    };
    /**
     * Computation error. Optional. Appears only if something happens on jobs execution.
     */
    errors?: {
      /**
       * Domain of the error appearance (internal).
       */
      domain?: "rupy" | "alfred" | "celim" | "webapp";
      /**
       * Should be a short, unique, machine-readable error code string. e.g. FileNotFound
       */
      reason?: string;
      /**
       * Human-readable error message. e.g. 'File Not Found: /home/demo/data/project1/job-123/job-config.json'
       */
      message?: string;
      /**
       * Full machine-readable error traceback. e.g. FileNotFound
       */
      traceback?: string;
    }[];
    /**
     * A Python compatible regex to exclude files from upload. e.g. ^.*.txt& excludes all files with .txt suffix
     */
    excludeFilesPattern?: string;
  };
  _project: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
  _material?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
  parent?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
  /**
   * Context variables that the job will have access to at runtime
   */
  runtimeContext?: {};
  /**
   * history of the workflow scope on each update
   */
  scopeTrack?: {
    repetition?: number;
    scope?: {
      global: {
        [k: string]: unknown;
      };
      local: {
        [k: string]: unknown;
      };
    };
  }[];
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  metadata?: {};
}
 
/** Schema dist/js/schema/job/compute.json */

/**
 * Custom keywords prefixed with validate correspond to custom validation methods implemented downstream
 */
export interface ComputeArgumentsSchema {
  /**
   * Name of the submission queues: https://docs.mat3ra.com/infrastructure/resource/queues/. Below enums are for Azure, then AWS circa 2022-08, hence the duplication.
   */
  queue:
    | "D"
    | "OR"
    | "OF"
    | "OFplus"
    | "SR"
    | "SF"
    | "SFplus"
    | "GPOF"
    | "GP2OF"
    | "GP4OF"
    | "GPSF"
    | "GP2SF"
    | "GP4SF"
    | "OR4"
    | "OR8"
    | "OR16"
    | "SR4"
    | "SR8"
    | "SR16"
    | "GOF"
    | "G4OF"
    | "G8OF"
    | "GSF"
    | "G4SF"
    | "G8SF";
  /**
   * number of nodes used for the job inside the RMS.
   */
  nodes: number;
  /**
   * number of CPUs used for the job inside the RMS.
   */
  ppn: number;
  /**
   * Wallclock time limit for computing a job. Clock format: 'hh:mm:ss'
   */
  timeLimit: string;
  /**
   * Convention to use when reasoning about time limits
   */
  timeLimitType?: "per single attempt" | "compound";
  /**
   * Job is allowed to restart on termination.
   */
  isRestartable?: boolean;
  /**
   * Email notification for the job: n - never, a - job aborted, b - job begins, e - job ends. Last three could be combined.
   */
  notify?: string;
  /**
   * Email address to notify about job execution.
   */
  email?: string;
  /**
   * Maximum CPU count per node. This parameter is used to let backend job submission infrastructure know that this job is to be charged for the maximum CPU per node instead of the actual ppn. For premium/fast queues where resources are provisioned on-demand and exclusively per user.
   */
  maxCPU?: number;
  /**
   * Optional arguments specific to using application - VASP, Quantum Espresso, etc. Specified elsewhere
   */
  arguments?: {
    /**
     * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
     */
    nimage?: number;
    /**
     * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
     */
    npools?: number;
    /**
     * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
     */
    nband?: number;
    /**
     * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
     */
    ntg?: number;
    /**
     * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
     */
    ndiag?: number;
  };
  /**
   * Cluster where the job is executed. Optional on create. Required on job submission.
   */
  cluster?: {
    /**
     * FQDN of the cluster. e.g. master-1-staging.exabyte.io
     */
    fqdn?: string;
    /**
     * Job's identity in RMS. e.g. 1234.master-1-staging.exabyte.io
     */
    jid?: string;
  };
  /**
   * Computation error. Optional. Appears only if something happens on jobs execution.
   */
  errors?: {
    /**
     * Domain of the error appearance (internal).
     */
    domain?: "rupy" | "alfred" | "celim" | "webapp";
    /**
     * Should be a short, unique, machine-readable error code string. e.g. FileNotFound
     */
    reason?: string;
    /**
     * Human-readable error message. e.g. 'File Not Found: /home/demo/data/project1/job-123/job-config.json'
     */
    message?: string;
    /**
     * Full machine-readable error traceback. e.g. FileNotFound
     */
    traceback?: string;
  }[];
  /**
   * A Python compatible regex to exclude files from upload. e.g. ^.*.txt& excludes all files with .txt suffix
   */
  excludeFilesPattern?: string;
}
 
/** Schema dist/js/schema/job.json */

export interface JobSchema {
  workflow: {
    /**
     * Array of subworkflows. Subworkflow can be an instance of workflow to allow for nesting
     */
    subworkflows: {
      /**
       * Contains the Units of the subworkflow
       */
      units: (
        | {
            /**
             * type of the unit
             */
            type: "io";
            subtype: "input" | "output" | "dataFrame";
            source: "api" | "db" | "object_storage";
            input: (
              | {
                  /**
                   * rest API endpoint
                   */
                  endpoint: string;
                  /**
                   * rest API endpoint options
                   */
                  endpoint_options: {};
                  /**
                   * the name of the variable in local scope to save the data under
                   */
                  name?: string;
                  [k: string]: unknown;
                }
              | (
                  | {
                      /**
                       * IDs of item to retrieve from db
                       */
                      ids: string[];
                      [k: string]: unknown;
                    }
                  | {
                      /**
                       * db collection name
                       */
                      collection: string;
                      /**
                       * whether the result should be saved as draft
                       */
                      draft: boolean;
                      [k: string]: unknown;
                    }
                )
              | {
                  objectData: {
                    /**
                     * Object storage container for the file
                     */
                    CONTAINER?: string;
                    /**
                     * Name of the file inside the object storage bucket
                     */
                    NAME?: string;
                    /**
                     * Object storage provider
                     */
                    PROVIDER?: string;
                    /**
                     * Region for the object container specified in Container
                     */
                    REGION?: string;
                    /**
                     * Size of the file in bytes
                     */
                    SIZE?: number;
                    /**
                     * Unix timestamp showing when the file was last modified
                     */
                    TIMESTAMP?: string;
                  };
                  /**
                   * if a file with the same filename already exists, whether to overwrite the old file
                   */
                  overwrite?: boolean;
                  /**
                   * Relative path to the directory that contains the file.
                   */
                  pathname?: string;
                  /**
                   * Basename of the file
                   */
                  basename?: string;
                  /**
                   * What kind of file this is, e.g. image / text
                   */
                  filetype?: string;
                  [k: string]: unknown;
                }
            )[];
            /**
             * entity identity
             */
            _id?: string;
            isDraft?: boolean;
            /**
             * name of the unit. e.g. pw_scf
             */
            name?: string;
            /**
             * Status of the unit.
             */
            status?: "idle" | "active" | "warning" | "error" | "finished";
            /**
             * Whether this unit is the first one to be executed.
             */
            head?: boolean;
            /**
             * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
             */
            flowchartId: string;
            /**
             * Next unit's flowchartId. If empty, the current unit is the last.
             */
            next?: string;
            /**
             * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
             */
            enableRender?: boolean;
            context?: {};
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * entity tags
             */
            tags?: string[];
            statusTrack?: {
              trackedAt: number;
              status: string;
              repetition?: number;
            }[];
            [k: string]: unknown;
          }
        | {
            /**
             * type of the unit
             */
            type: "reduce";
            /**
             * corresponding map unit flowchart ID
             */
            mapFlowchartId: string;
            /**
             * input information for reduce unit
             */
            input: {
              /**
               * reduce operation, e.g. aggregate
               */
              operation: string;
              /**
               * arguments which are passed to reduce operation function
               */
              arguments: string[];
            }[];
            /**
             * entity identity
             */
            _id?: string;
            isDraft?: boolean;
            /**
             * name of the unit. e.g. pw_scf
             */
            name?: string;
            /**
             * Status of the unit.
             */
            status?: "idle" | "active" | "warning" | "error" | "finished";
            /**
             * Whether this unit is the first one to be executed.
             */
            head?: boolean;
            /**
             * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
             */
            flowchartId: string;
            /**
             * Next unit's flowchartId. If empty, the current unit is the last.
             */
            next?: string;
            /**
             * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
             */
            enableRender?: boolean;
            context?: {};
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * entity tags
             */
            tags?: string[];
            statusTrack?: {
              trackedAt: number;
              status: string;
              repetition?: number;
            }[];
            [k: string]: unknown;
          }
        | {
            /**
             * type of the unit
             */
            type: "condition";
            /**
             * Input information for condition.
             */
            input: {
              /**
               * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
               */
              scope: string;
              /**
               * Name of the input data. e.g. total_energy
               */
              name: string;
            }[];
            /**
             * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
             */
            statement: string;
            /**
             * Flowchart ID reference for `then` part of the condition.
             */
            then: string;
            /**
             * Flowchart ID reference for `else` part of the condition.
             */
            else: string;
            /**
             * Maximum occurrence of the condition, usable for loops.
             */
            maxOccurrences: number;
            /**
             * Throw exception on reaching to maximum occurence.
             */
            throwException?: boolean;
            /**
             * entity identity
             */
            _id?: string;
            isDraft?: boolean;
            /**
             * name of the unit. e.g. pw_scf
             */
            name?: string;
            /**
             * Status of the unit.
             */
            status?: "idle" | "active" | "warning" | "error" | "finished";
            /**
             * Whether this unit is the first one to be executed.
             */
            head?: boolean;
            /**
             * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
             */
            flowchartId: string;
            /**
             * Next unit's flowchartId. If empty, the current unit is the last.
             */
            next?: string;
            /**
             * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
             */
            enableRender?: boolean;
            context?: {};
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * entity tags
             */
            tags?: string[];
            statusTrack?: {
              trackedAt: number;
              status: string;
              repetition?: number;
            }[];
            [k: string]: unknown;
          }
        | {
            /**
             * type of the unit
             */
            type: "assertion";
            /**
             * The statement to be evaluated
             */
            statement: string;
            /**
             * The error message to be displayed if the assertion fails
             */
            errorMessage?: string;
            /**
             * entity identity
             */
            _id?: string;
            isDraft?: boolean;
            /**
             * name of the unit. e.g. pw_scf
             */
            name: string;
            /**
             * Status of the unit.
             */
            status?: "idle" | "active" | "warning" | "error" | "finished";
            /**
             * Whether this unit is the first one to be executed.
             */
            head?: boolean;
            /**
             * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
             */
            flowchartId: string;
            /**
             * Next unit's flowchartId. If empty, the current unit is the last.
             */
            next?: string;
            /**
             * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
             */
            enableRender?: boolean;
            context?: {};
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * entity tags
             */
            tags?: string[];
            statusTrack?: {
              trackedAt: number;
              status: string;
              repetition?: number;
            }[];
            [k: string]: unknown;
          }
        | {
            /**
             * type of the unit
             */
            type: "execution";
            application: {
              /**
               * The short name of the application. e.g. qe
               */
              shortName?: string;
              /**
               * Application's short description.
               */
              summary?: string;
              /**
               * Application version. e.g. 5.3.5
               */
              version?: string;
              /**
               * Application build. e.g. VTST
               */
              build?: string;
              /**
               * Whether advanced compute options are present
               */
              hasAdvancedComputeOptions?: boolean;
              /**
               * Whether licensing is present
               */
              isLicensed?: boolean;
              /**
               * entity identity
               */
              _id?: string;
              /**
               * entity slug
               */
              slug?: string;
              systemName?: string;
              consistencyChecks?: {
                /**
                 * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
                 */
                key: string;
                /**
                 * Name of the consistency check that is performed, which is listed in an enum.
                 */
                name: "default" | "atomsTooClose" | "atomsOverlap";
                /**
                 * Severity level of the problem, which is used in UI to differentiate.
                 */
                severity: "info" | "warning" | "error";
                /**
                 * Message generated by the consistency check describing the problem.
                 */
                message: string;
              }[];
              /**
               * entity's schema version. Used to distinct between different schemas.
               */
              schemaVersion?: string;
              /**
               * entity name
               */
              name?: string;
              /**
               * Identifies that entity is defaultable
               */
              isDefault?: boolean;
              [k: string]: unknown;
            };
            executable?: {
              /**
               * The name of the executable. e.g. pw.x
               */
              name: string;
              /**
               * _ids of the application this executable belongs to
               */
              applicationId?: string[];
              /**
               * Whether advanced compute options are present
               */
              hasAdvancedComputeOptions?: boolean;
              /**
               * entity identity
               */
              _id?: string;
              /**
               * entity slug
               */
              slug?: string;
              systemName?: string;
              consistencyChecks?: {
                /**
                 * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
                 */
                key: string;
                /**
                 * Name of the consistency check that is performed, which is listed in an enum.
                 */
                name: "default" | "atomsTooClose" | "atomsOverlap";
                /**
                 * Severity level of the problem, which is used in UI to differentiate.
                 */
                severity: "info" | "warning" | "error";
                /**
                 * Message generated by the consistency check describing the problem.
                 */
                message: string;
              }[];
              /**
               * entity's schema version. Used to distinct between different schemas.
               */
              schemaVersion?: string;
              /**
               * Identifies that entity is defaultable
               */
              isDefault?: boolean;
              /**
               * names of the pre-processors for this calculation
               */
              preProcessors?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
              /**
               * names of the post-processors for this calculation
               */
              postProcessors?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
              /**
               * names of the monitors for this calculation
               */
              monitors?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
              /**
               * names of the results for this calculation
               */
              results?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
            };
            flavor?: {
              /**
               * _id of the executable this flavor belongs to
               */
              executableId?: string;
              /**
               * name of the executable this flavor belongs to
               */
              executableName?: string;
              /**
               * name of the application this flavor belongs to
               */
              applicationName?: string;
              input?: {
                templateId?: string;
                templateName?: string;
                /**
                 * name of the resulting input file, if different than template name
                 */
                name?: string;
              }[];
              /**
               * entity identity
               */
              _id?: string;
              /**
               * entity slug
               */
              slug?: string;
              systemName?: string;
              consistencyChecks?: {
                /**
                 * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
                 */
                key: string;
                /**
                 * Name of the consistency check that is performed, which is listed in an enum.
                 */
                name: "default" | "atomsTooClose" | "atomsOverlap";
                /**
                 * Severity level of the problem, which is used in UI to differentiate.
                 */
                severity: "info" | "warning" | "error";
                /**
                 * Message generated by the consistency check describing the problem.
                 */
                message: string;
              }[];
              /**
               * entity's schema version. Used to distinct between different schemas.
               */
              schemaVersion?: string;
              /**
               * entity name
               */
              name?: string;
              /**
               * Identifies that entity is defaultable
               */
              isDefault?: boolean;
              /**
               * names of the pre-processors for this calculation
               */
              preProcessors?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
              /**
               * names of the post-processors for this calculation
               */
              postProcessors?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
              /**
               * names of the monitors for this calculation
               */
              monitors?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
              /**
               * names of the results for this calculation
               */
              results?: (
                | {
                    /**
                     * The name of this item. e.g. scf_accuracy
                     */
                    name: string;
                  }
                | string
              )[];
            };
            /**
             * unit input (type to be specified by the application's execution unit)
             */
            input: {
              [k: string]: unknown;
            };
            /**
             * entity identity
             */
            _id?: string;
            isDraft?: boolean;
            /**
             * name of the unit. e.g. pw_scf
             */
            name?: string;
            /**
             * Status of the unit.
             */
            status?: "idle" | "active" | "warning" | "error" | "finished";
            /**
             * Whether this unit is the first one to be executed.
             */
            head?: boolean;
            /**
             * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
             */
            flowchartId: string;
            /**
             * Next unit's flowchartId. If empty, the current unit is the last.
             */
            next?: string;
            /**
             * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
             */
            enableRender?: boolean;
            context?: {};
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * entity tags
             */
            tags?: string[];
            statusTrack?: {
              trackedAt: number;
              status: string;
              repetition?: number;
            }[];
            [k: string]: unknown;
          }
        | {
            /**
             * type of the unit
             */
            type: "assignment";
            /**
             * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
             */
            input?: {
              /**
               * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
               */
              scope: string;
              /**
               * Name of the input data. e.g. total_energy
               */
              name: string;
            }[];
            /**
             * Name of the global variable. e.g. 'x'
             */
            operand: string;
            /**
             * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
             */
            value: string | boolean | number;
            /**
             * entity identity
             */
            _id?: string;
            isDraft?: boolean;
            /**
             * name of the unit. e.g. pw_scf
             */
            name: string;
            /**
             * Status of the unit.
             */
            status?: "idle" | "active" | "warning" | "error" | "finished";
            /**
             * Whether this unit is the first one to be executed.
             */
            head?: boolean;
            /**
             * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
             */
            flowchartId: string;
            /**
             * Next unit's flowchartId. If empty, the current unit is the last.
             */
            next?: string;
            /**
             * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
             */
            enableRender?: boolean;
            context?: {};
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * entity tags
             */
            tags?: string[];
            statusTrack?: {
              trackedAt: number;
              status: string;
              repetition?: number;
            }[];
            scope?: string;
            [k: string]: unknown;
          }
        | {
            /**
             * type of the unit
             */
            type: "processing";
            /**
             * Contains information about the operation used.
             */
            operation: string;
            /**
             * Contains information about the specific type of the operation used.
             */
            operationType: string;
            /**
             * unit input (type to be specified by the child units)
             */
            inputData: {
              [k: string]: unknown;
            };
            /**
             * entity identity
             */
            _id?: string;
            isDraft?: boolean;
            /**
             * name of the unit. e.g. pw_scf
             */
            name?: string;
            /**
             * Status of the unit.
             */
            status?: "idle" | "active" | "warning" | "error" | "finished";
            /**
             * Whether this unit is the first one to be executed.
             */
            head?: boolean;
            /**
             * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
             */
            flowchartId: string;
            /**
             * Next unit's flowchartId. If empty, the current unit is the last.
             */
            next?: string;
            /**
             * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
             */
            enableRender?: boolean;
            context?: {};
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * entity tags
             */
            tags?: string[];
            statusTrack?: {
              trackedAt: number;
              status: string;
              repetition?: number;
            }[];
            [k: string]: unknown;
          }
      )[];
      model: {
        /**
         * general type of the model, eg. `dft`
         */
        type: string;
        /**
         * general subtype of the model, eg. `lda`
         */
        subtype: string;
        method: {
          /**
           * general type of this method, eg. `pseudopotential`
           */
          type: string;
          /**
           * general subtype of this method, eg. `ultra-soft`
           */
          subtype: string;
          /**
           * Object showing the actual possible precision based on theory and implementation
           */
          precision?: {};
          /**
           * additional data specific to method, eg. array of pseudopotentials
           */
          data?: {};
        };
        [k: string]: unknown;
      };
      application: {
        /**
         * The short name of the application. e.g. qe
         */
        shortName?: string;
        /**
         * Application's short description.
         */
        summary?: string;
        /**
         * Application version. e.g. 5.3.5
         */
        version?: string;
        /**
         * Application build. e.g. VTST
         */
        build?: string;
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * Whether licensing is present
         */
        isLicensed?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        [k: string]: unknown;
      };
      /**
       * Defines whether to store the results/properties extracted in this unit to properties collection
       */
      isDraft?: boolean;
      /**
       * subworkflow identity
       */
      _id?: string;
      /**
       * Human-readable name of the subworkflow. e.g. Total-energy
       */
      name: string;
      /**
       * Array of characteristic properties calculated by this subworkflow
       */
      properties?: (string | {})[];
      /**
       * compute parameters
       */
      compute?: {
        /**
         * Name of the submission queues: https://docs.mat3ra.com/infrastructure/resource/queues/. Below enums are for Azure, then AWS circa 2022-08, hence the duplication.
         */
        queue:
          | "D"
          | "OR"
          | "OF"
          | "OFplus"
          | "SR"
          | "SF"
          | "SFplus"
          | "GPOF"
          | "GP2OF"
          | "GP4OF"
          | "GPSF"
          | "GP2SF"
          | "GP4SF"
          | "OR4"
          | "OR8"
          | "OR16"
          | "SR4"
          | "SR8"
          | "SR16"
          | "GOF"
          | "G4OF"
          | "G8OF"
          | "GSF"
          | "G4SF"
          | "G8SF";
        /**
         * number of nodes used for the job inside the RMS.
         */
        nodes: number;
        /**
         * number of CPUs used for the job inside the RMS.
         */
        ppn: number;
        /**
         * Wallclock time limit for computing a job. Clock format: 'hh:mm:ss'
         */
        timeLimit: string;
        /**
         * Convention to use when reasoning about time limits
         */
        timeLimitType?: "per single attempt" | "compound";
        /**
         * Job is allowed to restart on termination.
         */
        isRestartable?: boolean;
        /**
         * Email notification for the job: n - never, a - job aborted, b - job begins, e - job ends. Last three could be combined.
         */
        notify?: string;
        /**
         * Email address to notify about job execution.
         */
        email?: string;
        /**
         * Maximum CPU count per node. This parameter is used to let backend job submission infrastructure know that this job is to be charged for the maximum CPU per node instead of the actual ppn. For premium/fast queues where resources are provisioned on-demand and exclusively per user.
         */
        maxCPU?: number;
        /**
         * Optional arguments specific to using application - VASP, Quantum Espresso, etc. Specified elsewhere
         */
        arguments?: {
          /**
           * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
           */
          nimage?: number;
          /**
           * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
           */
          npools?: number;
          /**
           * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
           */
          nband?: number;
          /**
           * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
           */
          ntg?: number;
          /**
           * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
           */
          ndiag?: number;
        };
        /**
         * Cluster where the job is executed. Optional on create. Required on job submission.
         */
        cluster?: {
          /**
           * FQDN of the cluster. e.g. master-1-staging.exabyte.io
           */
          fqdn?: string;
          /**
           * Job's identity in RMS. e.g. 1234.master-1-staging.exabyte.io
           */
          jid?: string;
        };
        /**
         * Computation error. Optional. Appears only if something happens on jobs execution.
         */
        errors?: {
          /**
           * Domain of the error appearance (internal).
           */
          domain?: "rupy" | "alfred" | "celim" | "webapp";
          /**
           * Should be a short, unique, machine-readable error code string. e.g. FileNotFound
           */
          reason?: string;
          /**
           * Human-readable error message. e.g. 'File Not Found: /home/demo/data/project1/job-123/job-config.json'
           */
          message?: string;
          /**
           * Full machine-readable error traceback. e.g. FileNotFound
           */
          traceback?: string;
        }[];
        /**
         * A Python compatible regex to exclude files from upload. e.g. ^.*.txt& excludes all files with .txt suffix
         */
        excludeFilesPattern?: string;
      } | null;
    }[];
    /**
     * Contains the Units of the Workflow
     */
    units: (
      | {
          /**
           * type of the unit
           */
          type: "io";
          subtype: "input" | "output" | "dataFrame";
          source: "api" | "db" | "object_storage";
          input: (
            | {
                /**
                 * rest API endpoint
                 */
                endpoint: string;
                /**
                 * rest API endpoint options
                 */
                endpoint_options: {};
                /**
                 * the name of the variable in local scope to save the data under
                 */
                name?: string;
                [k: string]: unknown;
              }
            | (
                | {
                    /**
                     * IDs of item to retrieve from db
                     */
                    ids: string[];
                    [k: string]: unknown;
                  }
                | {
                    /**
                     * db collection name
                     */
                    collection: string;
                    /**
                     * whether the result should be saved as draft
                     */
                    draft: boolean;
                    [k: string]: unknown;
                  }
              )
            | {
                objectData: {
                  /**
                   * Object storage container for the file
                   */
                  CONTAINER?: string;
                  /**
                   * Name of the file inside the object storage bucket
                   */
                  NAME?: string;
                  /**
                   * Object storage provider
                   */
                  PROVIDER?: string;
                  /**
                   * Region for the object container specified in Container
                   */
                  REGION?: string;
                  /**
                   * Size of the file in bytes
                   */
                  SIZE?: number;
                  /**
                   * Unix timestamp showing when the file was last modified
                   */
                  TIMESTAMP?: string;
                };
                /**
                 * if a file with the same filename already exists, whether to overwrite the old file
                 */
                overwrite?: boolean;
                /**
                 * Relative path to the directory that contains the file.
                 */
                pathname?: string;
                /**
                 * Basename of the file
                 */
                basename?: string;
                /**
                 * What kind of file this is, e.g. image / text
                 */
                filetype?: string;
                [k: string]: unknown;
              }
          )[];
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "reduce";
          /**
           * corresponding map unit flowchart ID
           */
          mapFlowchartId: string;
          /**
           * input information for reduce unit
           */
          input: {
            /**
             * reduce operation, e.g. aggregate
             */
            operation: string;
            /**
             * arguments which are passed to reduce operation function
             */
            arguments: string[];
          }[];
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "condition";
          /**
           * Input information for condition.
           */
          input: {
            /**
             * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
             */
            scope: string;
            /**
             * Name of the input data. e.g. total_energy
             */
            name: string;
          }[];
          /**
           * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
           */
          statement: string;
          /**
           * Flowchart ID reference for `then` part of the condition.
           */
          then: string;
          /**
           * Flowchart ID reference for `else` part of the condition.
           */
          else: string;
          /**
           * Maximum occurrence of the condition, usable for loops.
           */
          maxOccurrences: number;
          /**
           * Throw exception on reaching to maximum occurence.
           */
          throwException?: boolean;
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "assertion";
          /**
           * The statement to be evaluated
           */
          statement: string;
          /**
           * The error message to be displayed if the assertion fails
           */
          errorMessage?: string;
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "execution";
          application: {
            /**
             * The short name of the application. e.g. qe
             */
            shortName?: string;
            /**
             * Application's short description.
             */
            summary?: string;
            /**
             * Application version. e.g. 5.3.5
             */
            version?: string;
            /**
             * Application build. e.g. VTST
             */
            build?: string;
            /**
             * Whether advanced compute options are present
             */
            hasAdvancedComputeOptions?: boolean;
            /**
             * Whether licensing is present
             */
            isLicensed?: boolean;
            /**
             * entity identity
             */
            _id?: string;
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * entity name
             */
            name?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            [k: string]: unknown;
          };
          executable?: {
            /**
             * The name of the executable. e.g. pw.x
             */
            name: string;
            /**
             * _ids of the application this executable belongs to
             */
            applicationId?: string[];
            /**
             * Whether advanced compute options are present
             */
            hasAdvancedComputeOptions?: boolean;
            /**
             * entity identity
             */
            _id?: string;
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
          };
          flavor?: {
            /**
             * _id of the executable this flavor belongs to
             */
            executableId?: string;
            /**
             * name of the executable this flavor belongs to
             */
            executableName?: string;
            /**
             * name of the application this flavor belongs to
             */
            applicationName?: string;
            input?: {
              templateId?: string;
              templateName?: string;
              /**
               * name of the resulting input file, if different than template name
               */
              name?: string;
            }[];
            /**
             * entity identity
             */
            _id?: string;
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * entity name
             */
            name?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
          };
          /**
           * unit input (type to be specified by the application's execution unit)
           */
          input: {
            [k: string]: unknown;
          };
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "assignment";
          /**
           * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
           */
          input?: {
            /**
             * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
             */
            scope: string;
            /**
             * Name of the input data. e.g. total_energy
             */
            name: string;
          }[];
          /**
           * Name of the global variable. e.g. 'x'
           */
          operand: string;
          /**
           * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
           */
          value: string | boolean | number;
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          scope?: string;
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "processing";
          /**
           * Contains information about the operation used.
           */
          operation: string;
          /**
           * Contains information about the specific type of the operation used.
           */
          operationType: string;
          /**
           * unit input (type to be specified by the child units)
           */
          inputData: {
            [k: string]: unknown;
          };
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "map";
          /**
           * Id of workflow to run inside map
           */
          workflowId: string;
          /**
           * Input information for map.
           */
          input: {
            /**
             * Name of the target variable to substitute using the values below. e.g. K_POINTS
             */
            target: string;
            /**
             * Scope to retrieve `values` from, global or flowchartId. Optional if `values` is given.
             */
            scope?: string;
            /**
             * Name of the variable inside the scope to retrieve `values` from. Optional if `values` is given.
             */
            name?: string;
            /**
             * Sequence of values for the target Jinja variable. Optional if `scope` and `name` are given. This can be used for map-reduce type parallel execution
             */
            values?: (string | number | {})[];
            useValues?: boolean;
          };
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "subworkflow";
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
    )[];
    /**
     * Array of characteristic properties calculated by this workflow (TODO: add enums)
     */
    properties?: (string | {})[];
    /**
     * Whether to use the dataset tab in the job designer. Mutually exclusive with using the materials tab.
     */
    isUsingDataset?: boolean;
    /**
     * Array of workflows with the same schema as the current one.
     */
    workflows?: {}[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    metadata?: {};
  };
  /**
   * Identity used to track jobs originated from command-line
   */
  rmsId?: string;
  /**
   * job status
   */
  status:
    | "pre-submission"
    | "queued"
    | "submitted"
    | "active"
    | "finished"
    | "terminate-queued"
    | "terminated"
    | "error"
    | "deleted"
    | "timeout";
  /**
   * Approximate start time of the job. e.g. within 10 min
   */
  startTime?: string;
  /**
   * The path to the working directory of this job, when the job originates from command-line
   */
  workDir?: string;
  /**
   * Custom keywords prefixed with validate correspond to custom validation methods implemented downstream
   */
  compute: {
    /**
     * Name of the submission queues: https://docs.mat3ra.com/infrastructure/resource/queues/. Below enums are for Azure, then AWS circa 2022-08, hence the duplication.
     */
    queue:
      | "D"
      | "OR"
      | "OF"
      | "OFplus"
      | "SR"
      | "SF"
      | "SFplus"
      | "GPOF"
      | "GP2OF"
      | "GP4OF"
      | "GPSF"
      | "GP2SF"
      | "GP4SF"
      | "OR4"
      | "OR8"
      | "OR16"
      | "SR4"
      | "SR8"
      | "SR16"
      | "GOF"
      | "G4OF"
      | "G8OF"
      | "GSF"
      | "G4SF"
      | "G8SF";
    /**
     * number of nodes used for the job inside the RMS.
     */
    nodes: number;
    /**
     * number of CPUs used for the job inside the RMS.
     */
    ppn: number;
    /**
     * Wallclock time limit for computing a job. Clock format: 'hh:mm:ss'
     */
    timeLimit: string;
    /**
     * Convention to use when reasoning about time limits
     */
    timeLimitType?: "per single attempt" | "compound";
    /**
     * Job is allowed to restart on termination.
     */
    isRestartable?: boolean;
    /**
     * Email notification for the job: n - never, a - job aborted, b - job begins, e - job ends. Last three could be combined.
     */
    notify?: string;
    /**
     * Email address to notify about job execution.
     */
    email?: string;
    /**
     * Maximum CPU count per node. This parameter is used to let backend job submission infrastructure know that this job is to be charged for the maximum CPU per node instead of the actual ppn. For premium/fast queues where resources are provisioned on-demand and exclusively per user.
     */
    maxCPU?: number;
    /**
     * Optional arguments specific to using application - VASP, Quantum Espresso, etc. Specified elsewhere
     */
    arguments?: {
      /**
       * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
       */
      nimage?: number;
      /**
       * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
       */
      npools?: number;
      /**
       * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
       */
      nband?: number;
      /**
       * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
       */
      ntg?: number;
      /**
       * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
       */
      ndiag?: number;
    };
    /**
     * Cluster where the job is executed. Optional on create. Required on job submission.
     */
    cluster?: {
      /**
       * FQDN of the cluster. e.g. master-1-staging.exabyte.io
       */
      fqdn?: string;
      /**
       * Job's identity in RMS. e.g. 1234.master-1-staging.exabyte.io
       */
      jid?: string;
    };
    /**
     * Computation error. Optional. Appears only if something happens on jobs execution.
     */
    errors?: {
      /**
       * Domain of the error appearance (internal).
       */
      domain?: "rupy" | "alfred" | "celim" | "webapp";
      /**
       * Should be a short, unique, machine-readable error code string. e.g. FileNotFound
       */
      reason?: string;
      /**
       * Human-readable error message. e.g. 'File Not Found: /home/demo/data/project1/job-123/job-config.json'
       */
      message?: string;
      /**
       * Full machine-readable error traceback. e.g. FileNotFound
       */
      traceback?: string;
    }[];
    /**
     * A Python compatible regex to exclude files from upload. e.g. ^.*.txt& excludes all files with .txt suffix
     */
    excludeFilesPattern?: string;
  };
  _project: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
  _material?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
  parent?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
  /**
   * Context variables that the job will have access to at runtime
   */
  runtimeContext?: {};
  /**
   * history of the workflow scope on each update
   */
  scopeTrack?: {
    repetition?: number;
    scope?: {
      global: {
        [k: string]: unknown;
      };
      local: {
        [k: string]: unknown;
      };
    };
  }[];
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  metadata?: {};
}
 
/** Schema dist/js/schema/material/conventional.json */

export interface MaterialConventionalSchema {
  conventional?: {};
}
 
/** Schema dist/js/schema/material.json */

export interface MaterialSchema {
  /**
   * reduced chemical formula
   */
  formula?: string;
  /**
   * chemical formula based on the number of atoms of each element in the supercell
   */
  unitCellFormula?: string;
  basis: {
    elements: {
      id: number;
      value: string;
      /**
       * Occurrence is for fractional occupations
       */
      occurrence?: number;
      oxidationState?: number;
    }[];
    /**
     * Optional numeric label (e.g., 1, 2, as in Fe1, Fe2) to distinguish same atomic species to attach different spin magnetic moment.
     */
    labels?: {
      id?: number;
      value?: number;
    }[];
    coordinates: {
      id?: number;
      value?: [number, number, number] | [boolean, boolean, boolean];
    }[];
    name?: string;
    units?: string;
    bonds?: {
      /**
       * indices of the two connected atoms
       *
       * @minItems 2
       * @maxItems 2
       */
      atomPair?: [
        {
          /**
           * integer id of this entry
           */
          id?: number;
        },
        {
          /**
           * integer id of this entry
           */
          id?: number;
        }
      ];
      bondType?: "single" | "double" | "triple" | "quadruple" | "aromatic" | "tautomeric" | "dative" | "other";
    }[];
  };
  lattice: {
    name?: "lattice";
    vectors?: {
      /**
       * lattice parameter for fractional coordinates
       */
      alat?: number;
      units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
      /**
       * @minItems 3
       * @maxItems 3
       */
      a: [number, number, number];
      /**
       * @minItems 3
       * @maxItems 3
       */
      b: [number, number, number];
      /**
       * @minItems 3
       * @maxItems 3
       */
      c: [number, number, number];
    };
    type:
      | "CUB"
      | "BCC"
      | "FCC"
      | "TET"
      | "MCL"
      | "ORC"
      | "ORCC"
      | "ORCF"
      | "ORCI"
      | "HEX"
      | "BCT"
      | "TRI"
      | "MCLC"
      | "RHL";
    units?: {
      length?: "angstrom" | "bohr";
      angle?: "degree" | "radian";
    };
    /**
     * length of the first lattice vector
     */
    a: number;
    /**
     * length of the second lattice vector
     */
    b: number;
    /**
     * length of the third lattice vector
     */
    c: number;
    /**
     * angle between first and second lattice vector
     */
    alpha: number;
    /**
     * angle between second and third lattice vector
     */
    beta: number;
    /**
     * angle between first and third lattice vector
     */
    gamma: number;
  };
  derivedProperties?: (
    | {
        name?: "volume";
        units?: "angstrom^3";
        value: number;
      }
    | {
        name?: "density";
        units?: "g/cm^3";
        value: number;
      }
    | {
        /**
         * point group symbol in Schoenflies notation
         */
        pointGroupSymbol?: string;
        /**
         * space group symbol in Hermann–Mauguin notation
         */
        spaceGroupSymbol?: string;
        /**
         * tolerance used for symmetry calculation
         */
        tolerance?: {
          units?: "angstrom";
          value: number;
        };
        name?: "symmetry";
      }
    | {
        name?: "elemental_ratio";
        value: number;
        /**
         * the element this ratio is for
         */
        element?: string;
      }
    | {
        name?: "p-norm";
        /**
         * degree of the dimensionality of the norm
         */
        degree?: number;
        value: number;
      }
    | {
        name?: "inchi";
        value: string;
      }
    | {
        name?: "inchi_key";
        value: string;
      }
  )[];
  /**
   * information about a database source
   */
  external?: {
    /**
     * ID string for the materials uploaded from a third party source inside the third party source. For materialsproject.org an example ID is mp-32
     */
    id: string | number;
    /**
     * Third party source name, e.g. materials project, 2dmatpedia, ICSD, etc.
     */
    source: string;
    /**
     * Deprecated. To be removed. A flag that is true when material is initially imported from a third party * (as opposed to being independently designed from scratch).
     */
    origin: boolean;
    /**
     * Original response from external source.
     */
    data?: {};
    /**
     * Digital Object Identifier, e.g. 10.1088/0953-8984/25/10/105506
     */
    doi?: string;
    /**
     * The URL of the original record, e.g. https://next-gen.materialsproject.org/materials/mp-48; ToDo: update to use URI type per https://json-schema.org/understanding-json-schema/reference/string#resource-identifiers
     */
    url?: string;
  };
  /**
   * file source with the information inside
   */
  src?: {
    /**
     * file extension
     */
    extension?: string;
    /**
     * file name without extension
     */
    filename: string;
    /**
     * file content as raw text
     */
    text: string;
    /**
     * MD5 hash based on file content
     */
    hash: string;
  };
  /**
   * Hash string for a scaled structure with lattice vector a set to 1 (eg. for materials under pressure).
   */
  scaledHash?: string;
  /**
   * Corresponding ICSD id of the material
   */
  icsdId?: number;
  /**
   * Whether to work in the finite molecular picture (usually with atomic orbital basis)
   */
  isNonPeriodic?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  metadata?: {};
}
 
/** Schema dist/js/schema/method/categorized_method.json */

export interface CategorizedMethod {
  units: {
    /**
     * Used to categorize entities such as models and methods
     */
    categories?: {
      /**
       * contains either object with slugified entry or slug only as a string
       */
      tier1?:
        | {
            /**
             * descriptive human-readable name of entry
             */
            name: string;
            /**
             * machine-readable identifier
             */
            slug: string;
          }
        | string;
      /**
       * contains either object with slugified entry or slug only as a string
       */
      tier2?:
        | {
            /**
             * descriptive human-readable name of entry
             */
            name: string;
            /**
             * machine-readable identifier
             */
            slug: string;
          }
        | string;
      /**
       * contains either object with slugified entry or slug only as a string
       */
      tier3?:
        | {
            /**
             * descriptive human-readable name of entry
             */
            name: string;
            /**
             * machine-readable identifier
             */
            slug: string;
          }
        | string;
      /**
       * contains either object with slugified entry or slug only as a string
       */
      type?:
        | {
            /**
             * descriptive human-readable name of entry
             */
            name: string;
            /**
             * machine-readable identifier
             */
            slug: string;
          }
        | string;
      /**
       * contains either object with slugified entry or slug only as a string
       */
      subtype?:
        | {
            /**
             * descriptive human-readable name of entry
             */
            name: string;
            /**
             * machine-readable identifier
             */
            slug: string;
          }
        | string;
    };
    /**
     * Instructive parameters defining the method
     */
    parameters?: {};
    /**
     * Object showing the actual possible precision based on theory and implementation
     */
    precision?: {};
    /**
     * entity name
     */
    name?: string;
    /**
     * TODO: Use regex once schema draft version has been updated
     */
    path?: string;
    /**
     * entity tags
     */
    tags?: string[];
  }[];
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/method/method_parameters.json */

export type MethodParameters =
  | {
      basisSlug?: "cc-pvdz" | "cc-pvtz" | "cc-pvqz";
    }
  | {
      basisSlug?: "3-21G" | "6-31G" | "6-311G";
    }
  | {
      basisSlug?: "sto-3g" | "sto-4g" | "sto-6g" | "def2-svp" | "def2-tzvp" | "def2-qzvp" | "cbs-qb3";
    };
 
/** Schema dist/js/schema/method/unit_method.json */

export interface CategorizedUnitMethod {
  /**
   * Used to categorize entities such as models and methods
   */
  categories?: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/method.json */

export interface BaseMethod {
  /**
   * general type of this method, eg. `pseudopotential`
   */
  type: string;
  /**
   * general subtype of this method, eg. `ultra-soft`
   */
  subtype: string;
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * additional data specific to method, eg. array of pseudopotentials
   */
  data?: {};
}
 
/** Schema dist/js/schema/methods_category/mathematical/diff/enum_options.json */

export interface MethodsCategoryMathematicalDiffEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/diff/fd.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface FiniteDifferenceMethodCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fd";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/diff.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface NumericalDifferentiationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/enum_options.json */

export interface MethodsCategoryMathematicalDiscrEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/mesh/enum_options.json */

export interface MethodsCategoryMathematicalDiscrMeshEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/mesh/hybrid.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface HybridMeshingCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "hybrid";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "mesh";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "discr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/mesh/nstruct.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface UnstructuredMeshingCategoryNstructSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "nstruct";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "mesh";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "discr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/mesh/struct/cartesian.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface CartesianGridSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "cartesian";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "struct";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "mesh";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "discr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/mesh/struct/enum_options.json */

export interface MethodsCategoryMathematicalDiscrMeshStructEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/mesh/struct.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface StructuredMeshingCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "struct";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "mesh";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "discr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr/mesh.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MeshingMethodCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "mesh";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "discr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/discr.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DiscretizationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "discr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/enum_options.json */

export interface MethodsCategoryMathematicalEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx/basisexp.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface BasisExpansionCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "basisExp";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fapprx";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx/enum_options.json */

export interface MethodsCategoryMathematicalFapprxEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx/ipol/enum_options.json */

export interface MethodsCategoryMathematicalFapprxIpolEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx/ipol/lin.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface LinearInterpolationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "lin";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ipol";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fapprx";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx/ipol/poly.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface PolynomialInterpolationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "poly";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ipol";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fapprx";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx/ipol/spline.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface SplineInterpolationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "spline";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ipol";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fapprx";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx/ipol.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface InterpolationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ipol";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fapprx";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/fapprx.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface UnstructuredMeshingCategoryFapprxSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fapprx";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/analytic/enum_options.json */

export interface MethodsCategoryMathematicalIntgrAnalyticEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/analytic/volume.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface AnalyticVolumeIntegralCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "volume";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("sphere" | "cube" | "rect-prism" | "tri-prism" | "cylinder" | "cone" | "tetrahedron" | "sq-pyr");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "analytic";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/analytic.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface AnalyticIntegralCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "analytic";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/diffeq/enum_options.json */

export interface MethodsCategoryMathematicalIntgrDiffeqEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/diffeq/order1.json */

/**
 * Categories for the numerical integration of differential equations
 */
export interface Order1Schema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "order1";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diffeq";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/diffeq/order2.json */

/**
 * Categories for the numerical integration of differential equations
 */
export interface Order2Schema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "order2";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diffeq";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/diffeq.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MethodsForTheNumericalIntegrationOfDifferentialEquationsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diffeq";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/enum_options.json */

export interface MethodsCategoryMathematicalIntgrEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/numquad/enum_options.json */

export interface MethodsCategoryMathematicalIntgrNumquadEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/numquad/gauss.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface GaussianQuadratureRulesSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "gauss";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "numquad";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/numquad/newcot.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface NewtonCotesQuadratureRulesSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "newcot";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "numquad";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/numquad.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MethodsForTheNumericalQuadratureSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "numquad";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/transf/enum_options.json */

export interface MethodsCategoryMathematicalIntgrTransfEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/transf/fourier.json */

/**
 * Fourier transform methods
 */
export interface FourierTransformMethodsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "fourier";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "transf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr/transf.json */

/**
 * Integral transform methods
 */
export interface IntegralTransformMethodsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "transf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/intgr.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface IntegrationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "intgr";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg/dcomp.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MatrixDecompositionMethodsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dcomp";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "linalg";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg/diag/davidson.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DavidsonDiagonalizationMethodSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "davidson";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diag";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "linalg";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg/diag/enum_options.json */

export interface MethodsCategoryMathematicalLinalgDiagEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg/diag.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MatrixDiagonalizationMethodsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diag";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "linalg";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg/enum_options.json */

export interface MethodsCategoryMathematicalLinalgEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg/lintra.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface LinearTransformationMethodsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "lintra";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "linalg";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg/matf.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MatrixFunctionMethodsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "matf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "linalg";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/linalg.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface LinearAlgebraCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "linalg";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/bracket.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface BracketAlgorithmsForTheOptimizationOfDifferentiableFunctionsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "bracket";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/enum_options.json */

export interface MethodsCategoryMathematicalOptDiffEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/local.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface LocalDescentMethodsForTheOptimizationOfDifferentiableFunctionsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "local";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/order1.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface FirstOrderAlgorithmsForTheOptimizationOfDifferentiableFunctionsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "order1";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/order2.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface SecondOrderAlgorithmsForTheOptimizationOfDifferentiableFunctionsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "order2";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/ordern/cg.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface ConjugateGradientMethodSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "cg";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ordern";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/ordern/enum_options.json */

export interface MethodsCategoryMathematicalOptDiffOrdernEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff/ordern.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MixedOrderAndHigherOrderAlgorithmsForTheOptimizationOfDifferentiableFunctionsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ordern";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/diff.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface OptimizationMethodsForDifferentiableFunctionsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "diff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/enum_options.json */

export interface MethodsCategoryMathematicalOptEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/ndiff/direct.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DirectAlgorithmsForTheOptimizationOfNonDifferentiableFunctionsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "direct";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ndiff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/ndiff/enum_options.json */

export interface MethodsCategoryMathematicalOptNdiffEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/ndiff/pop.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface PopulationAlgorithmsForTheOptmizationOfNonDifferentiableFunctionsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pop";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ndiff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/ndiff/stoch.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface StochasticAlgorithmsForTheOptmizationOfNonDifferentiableFunctionsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "stoch";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ndiff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/ndiff.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface OptimizationMethodsForNonDifferentiableFunctionsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ndiff";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/root/bracket.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface BracketingMethodForFindingRootsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "bracket";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "root";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/root/enum_options.json */

export interface MethodsCategoryMathematicalOptRootEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/root/iter.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface IterativeMethodForRootFindingCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "iterative";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "root";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt/root.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface RootFindingCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "root";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/opt.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MathematicalOptSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "opt";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/mathematical/regression.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface LinearMethodsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("linear" | "kernel_ridge");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("least_squares" | "ridge");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/enum_options.json */

export interface MethodsCategoryPhysicalEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/enum_options.json */

export interface MethodsCategoryPhysicalQmEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/ao/dunning.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DunningCorrelationConsistentBasisSetCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dunning";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ao";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/ao/other.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface OtherNeitherPopleNorDunningBasisSetCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "other";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ao";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/ao/pople.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface PopleBasisSetCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pople";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ao";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/ao.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface ApproximatingTheElectronicWaveFunctionWithAAtomicOrbitalBasisSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ao";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("pople" | "dunning" | "other");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/enum_options.json */

export interface MethodsCategoryPhysicalQmWfEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/psp.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface PseudopotentialCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "psp";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("us" | "nc" | "paw" | "coulomb");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/pw.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface PlaneWaveCatgeorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pw";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/smearing.json */

/**
 * Approximating Heaviside step function with smooth function
 */
export interface SmearingMethodsCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "smearing";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("gaussian" | "marzari-vanderbilt" | "methfessel-paxton" | "fermi-dirac");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf/tetrahedron.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface TetrahedronMethodForBrillouinZoneIntegrationCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "tetrahedron";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("linear" | "optimized" | "bloechl");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm/wf.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MethodsRelatedToWaveFunctionsSchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "wf";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_category/physical/qm.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface QuantumMechanicalMethodCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/methods_directory/legacy/localorbital.json */

export interface LegacyMethodLocalorbital {
  /**
   * general type of this method, eg. `pseudopotential`
   */
  type: "localorbital";
  /**
   * general subtype of this method, eg. `ultra-soft`
   */
  subtype: "pople";
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * additional data specific to method, eg. array of pseudopotentials
   */
  data?: {};
}
 
/** Schema dist/js/schema/methods_directory/legacy/pseudopotential.json */

export interface LegacyMethodPseudopotential {
  /**
   * general type of this method, eg. `pseudopotential`
   */
  type: "pseudopotential";
  /**
   * general subtype of this method, eg. `ultra-soft`
   */
  subtype: "paw" | "nc" | "us" | "any";
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * additional data specific to method, eg. array of pseudopotentials
   */
  data?: {};
}
 
/** Schema dist/js/schema/methods_directory/legacy/regression.json */

export interface LegacyMethodRegression {
  /**
   * general type of this method, eg. `pseudopotential`
   */
  type: "linear" | "kernel_ridge";
  /**
   * general subtype of this method, eg. `ultra-soft`
   */
  subtype: "least_squares" | "ridge";
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision: {
    perProperty?: {
      /**
       * property name in 'flattened' format
       */
      name?: string;
      /**
       * training error of the estimator
       */
      trainingError: number;
      /**
       * prediction score of the estimator. Eg: r2_score
       */
      score?: number;
    }[];
  };
  /**
   * additional data specific to method, eg. array of pseudopotentials
   */
  data: {
    perProperty?: (
      | {
          /**
           * intercept (shift) from the linear or non-linear fit of data points
           */
          intercept: number;
          /**
           * per-feature (property used for training the ML method/model) parameters
           */
          perFeature: {
            /**
             * coefficient in linear regression
             */
            coefficient?: number;
            /**
             * feature name
             */
            name: string;
            /**
             * pvalue: https://en.wikipedia.org/wiki/P-value
             */
            importance?: number;
          }[];
        }
      | {
          /**
           * training data
           */
          xFit: unknown[];
          /**
           * dual coefficients
           */
          dualCoefficients: unknown[];
          /**
           * per-feature (property used for training the ML method/model) parameters
           */
          perFeature: {
            /**
             * coefficient in linear regression
             */
            coefficient?: number;
            /**
             * feature name
             */
            name: string;
            /**
             * pvalue: https://en.wikipedia.org/wiki/P-value
             */
            importance?: number;
          }[];
        }
    )[];
    /**
     * dataset for ml
     */
    dataSet?: {
      /**
       * array of exabyteIds for materials in dataset
       */
      exabyteIds: string[];
      /**
       * holder for any extra information, eg. coming from user-uploaded CSV file
       */
      extra?: {
        [k: string]: unknown;
      };
    };
  };
}
 
/** Schema dist/js/schema/methods_directory/legacy/unknown.json */

export interface LegacyMethodUnknown {
  /**
   * general type of this method, eg. `pseudopotential`
   */
  type: "unknown";
  /**
   * general subtype of this method, eg. `ultra-soft`
   */
  subtype: "unknown";
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * additional data specific to method, eg. array of pseudopotentials
   */
  data?: {};
}
 
/** Schema dist/js/schema/methods_directory/mathematical/cg.json */

/**
 * conjugate gradient method schema
 */
export interface UnitMethodConjugateGradient {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "cg";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ordern";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "diff";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "opt";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/methods_directory/mathematical/davidson.json */

/**
 * Davidson diagonalization method
 */
export interface UnitMethodDavidsonSchema {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "davidson";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "diag";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "linalg";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression/data.json */

export interface RegressionData {
  perProperty?: (
    | {
        /**
         * intercept (shift) from the linear or non-linear fit of data points
         */
        intercept: number;
        /**
         * per-feature (property used for training the ML method/model) parameters
         */
        perFeature: {
          /**
           * coefficient in linear regression
           */
          coefficient?: number;
          /**
           * feature name
           */
          name: string;
          /**
           * pvalue: https://en.wikipedia.org/wiki/P-value
           */
          importance?: number;
        }[];
      }
    | {
        /**
         * training data
         */
        xFit: unknown[];
        /**
         * dual coefficients
         */
        dualCoefficients: unknown[];
        /**
         * per-feature (property used for training the ML method/model) parameters
         */
        perFeature: {
          /**
           * coefficient in linear regression
           */
          coefficient?: number;
          /**
           * feature name
           */
          name: string;
          /**
           * pvalue: https://en.wikipedia.org/wiki/P-value
           */
          importance?: number;
        }[];
      }
  )[];
  /**
   * dataset for ml
   */
  dataSet?: {
    /**
     * array of exabyteIds for materials in dataset
     */
    exabyteIds: string[];
    /**
     * holder for any extra information, eg. coming from user-uploaded CSV file
     */
    extra?: {
      [k: string]: unknown;
    };
  };
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression/dataset.json */

/**
 * dataset for ml
 */
export interface MethodsDirectoryMathematicalRegressionDataset {
  /**
   * array of exabyteIds for materials in dataset
   */
  exabyteIds: string[];
  /**
   * holder for any extra information, eg. coming from user-uploaded CSV file
   */
  extra?: {
    [k: string]: unknown;
  };
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression/kernel_ridge/data_per_property.json */

export interface KernelRidgeRegressionParametersSchema {
  /**
   * training data
   */
  xFit: unknown[];
  /**
   * dual coefficients
   */
  dualCoefficients: unknown[];
  /**
   * per-feature (property used for training the ML method/model) parameters
   */
  perFeature: {
    /**
     * coefficient in linear regression
     */
    coefficient?: number;
    /**
     * feature name
     */
    name: string;
    /**
     * pvalue: https://en.wikipedia.org/wiki/P-value
     */
    importance?: number;
  }[];
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression/linear/data_per_property.json */

export interface LinearRegressionParametersSchema {
  /**
   * intercept (shift) from the linear or non-linear fit of data points
   */
  intercept: number;
  /**
   * per-feature (property used for training the ML method/model) parameters
   */
  perFeature: {
    /**
     * coefficient in linear regression
     */
    coefficient?: number;
    /**
     * feature name
     */
    name: string;
    /**
     * pvalue: https://en.wikipedia.org/wiki/P-value
     */
    importance?: number;
  }[];
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression/per_feature_item.json */

export interface PerFeaturePropertyUsedForTrainingTheMLMethodModelParametersSchema {
  /**
   * coefficient in linear regression
   */
  coefficient?: number;
  /**
   * feature name
   */
  name: string;
  /**
   * pvalue: https://en.wikipedia.org/wiki/P-value
   */
  importance?: number;
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression/precision.json */

export interface RegressionPrecision {
  perProperty?: {
    /**
     * property name in 'flattened' format
     */
    name?: string;
    /**
     * training error of the estimator
     */
    trainingError: number;
    /**
     * prediction score of the estimator. Eg: r2_score
     */
    score?: number;
  }[];
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression/precision_per_property.json */

export interface RegressionPrecisionPerPropertySchema {
  /**
   * property name in 'flattened' format
   */
  name?: string;
  /**
   * training error of the estimator
   */
  trainingError: number;
  /**
   * prediction score of the estimator. Eg: r2_score
   */
  score?: number;
}
 
/** Schema dist/js/schema/methods_directory/mathematical/regression.json */

export interface UnitMethodRegression {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      ("linear" | "kernel_ridge");
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      ("least_squares" | "ridge");
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision: {
    perProperty?: {
      /**
       * property name in 'flattened' format
       */
      name?: string;
      /**
       * training error of the estimator
       */
      trainingError: number;
      /**
       * prediction score of the estimator. Eg: r2_score
       */
      score?: number;
    }[];
  };
  data: {
    perProperty?: (
      | {
          /**
           * intercept (shift) from the linear or non-linear fit of data points
           */
          intercept: number;
          /**
           * per-feature (property used for training the ML method/model) parameters
           */
          perFeature: {
            /**
             * coefficient in linear regression
             */
            coefficient?: number;
            /**
             * feature name
             */
            name: string;
            /**
             * pvalue: https://en.wikipedia.org/wiki/P-value
             */
            importance?: number;
          }[];
        }
      | {
          /**
           * training data
           */
          xFit: unknown[];
          /**
           * dual coefficients
           */
          dualCoefficients: unknown[];
          /**
           * per-feature (property used for training the ML method/model) parameters
           */
          perFeature: {
            /**
             * coefficient in linear regression
             */
            coefficient?: number;
            /**
             * feature name
             */
            name: string;
            /**
             * pvalue: https://en.wikipedia.org/wiki/P-value
             */
            importance?: number;
          }[];
        }
    )[];
    /**
     * dataset for ml
     */
    dataSet?: {
      /**
       * array of exabyteIds for materials in dataset
       */
      exabyteIds: string[];
      /**
       * holder for any extra information, eg. coming from user-uploaded CSV file
       */
      extra?: {
        [k: string]: unknown;
      };
    };
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/methods_directory/physical/ao/dunning.json */

/**
 * Dunning correlation-consistent basis set unit method
 */
export interface UnitMethodAoDunning {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "dunning";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ao";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "wf";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {
    basisSlug?: "cc-pvdz" | "cc-pvtz" | "cc-pvqz";
  };
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
/**
 * This interface was referenced by `UnitMethodAoDunning`'s JSON-Schema
 * via the `definition` "ao-basis-dunning".
 */
export interface AoBasisDunning {
  basisSlug?: "cc-pvdz" | "cc-pvtz" | "cc-pvqz";
}
 
/** Schema dist/js/schema/methods_directory/physical/ao/enum_options.json */

export interface MethodsDirectoryPhysicalAoEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/methods_directory/physical/ao/other.json */

/**
 * Other (neither Pople nor Dunning) basis set unit method
 */
export interface UnitMethodAoOther {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "other";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ao";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "wf";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {
    basisSlug?: "sto-3g" | "sto-4g" | "sto-6g" | "def2-svp" | "def2-tzvp" | "def2-qzvp" | "cbs-qb3";
  };
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
/**
 * This interface was referenced by `UnitMethodAoOther`'s JSON-Schema
 * via the `definition` "ao-basis-other".
 */
export interface AoBasisOther {
  basisSlug?: "sto-3g" | "sto-4g" | "sto-6g" | "def2-svp" | "def2-tzvp" | "def2-qzvp" | "cbs-qb3";
}
 
/** Schema dist/js/schema/methods_directory/physical/ao/pople.json */

/**
 * Pople basis set unit method
 */
export interface UnitMethodAoPople {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pople";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ao";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "wf";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {
    basisSlug?: "3-21G" | "6-31G" | "6-311G";
  };
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
/**
 * This interface was referenced by `UnitMethodAoPople`'s JSON-Schema
 * via the `definition` "ao-basis-pople".
 */
export interface AoBasisPople {
  basisSlug?: "3-21G" | "6-31G" | "6-311G";
}
 
/** Schema dist/js/schema/methods_directory/physical/psp/file.json */

export interface PseudopotentialFile {
  slug?: "pseudopotential";
  data?: {
    /**
     * chemical element
     */
    element: string;
    /**
     * MD5 hash of the pseudopotential file
     */
    hash?: string;
    type: "us" | "nc" | "paw" | "coulomb";
    /**
     * explains where this came from
     */
    source: string;
    /**
     * explains the version of where this came from
     */
    version?: string;
    exchangeCorrelation: {
      /**
       * DFT approximation
       */
      approximation?: string;
      /**
       * Exchange correlation functional
       */
      functional?: string;
      /**
       * TODO: Use regex once schema draft version has been updated
       */
      path?: string;
    };
    /**
     * contains pseudo orbital information, including orbital names and occupations
     */
    valenceConfiguration?: {
      orbitalName?: string;
      orbitalIndex?: number;
      principalNumber?: number;
      angularMomentum?: number;
      /**
       * Shell occupation
       */
      occupation?: number;
    }[];
    /**
     * location of the pseudopotential file on filesystem
     */
    path: string;
    /**
     * The names of the simulation engines that can use this pseudopotential, e.g. espresso
     */
    apps: string[];
    /**
     * filename of pseudopotential file on filesystem
     */
    filename?: string;
    /**
     * name of the data category
     */
    name?: "pseudopotential";
  };
  /**
   * TODO: remove in the future
   */
  source?: {
    info?: {};
    type?: string;
  };
}
 
/** Schema dist/js/schema/methods_directory/physical/psp/file_data_item.json */

export interface FileDataItem {
  /**
   * chemical element
   */
  element: string;
  /**
   * MD5 hash of the pseudopotential file
   */
  hash?: string;
  type: "us" | "nc" | "paw" | "coulomb";
  /**
   * explains where this came from
   */
  source: string;
  /**
   * explains the version of where this came from
   */
  version?: string;
  exchangeCorrelation: {
    /**
     * DFT approximation
     */
    approximation?: string;
    /**
     * Exchange correlation functional
     */
    functional?: string;
    /**
     * TODO: Use regex once schema draft version has been updated
     */
    path?: string;
  };
  /**
   * contains pseudo orbital information, including orbital names and occupations
   */
  valenceConfiguration?: {
    orbitalName?: string;
    orbitalIndex?: number;
    principalNumber?: number;
    angularMomentum?: number;
    /**
     * Shell occupation
     */
    occupation?: number;
  }[];
  /**
   * location of the pseudopotential file on filesystem
   */
  path: string;
  /**
   * The names of the simulation engines that can use this pseudopotential, e.g. espresso
   */
  apps: string[];
  /**
   * filename of pseudopotential file on filesystem
   */
  filename?: string;
  /**
   * name of the data category
   */
  name?: "pseudopotential";
}
 
/** Schema dist/js/schema/methods_directory/physical/psp.json */

/**
 * Core-valence separation by means of pseudopotentials (effective potential)
 */
export interface UnitMethodPseudopotential {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "psp";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      ("us" | "nc" | "paw" | "coulomb");
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "wf";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  data?: {
    slug?: "pseudopotential";
    data?: {
      /**
       * chemical element
       */
      element: string;
      /**
       * MD5 hash of the pseudopotential file
       */
      hash?: string;
      type: "us" | "nc" | "paw" | "coulomb";
      /**
       * explains where this came from
       */
      source: string;
      /**
       * explains the version of where this came from
       */
      version?: string;
      exchangeCorrelation: {
        /**
         * DFT approximation
         */
        approximation?: string;
        /**
         * Exchange correlation functional
         */
        functional?: string;
        /**
         * TODO: Use regex once schema draft version has been updated
         */
        path?: string;
      };
      /**
       * contains pseudo orbital information, including orbital names and occupations
       */
      valenceConfiguration?: {
        orbitalName?: string;
        orbitalIndex?: number;
        principalNumber?: number;
        angularMomentum?: number;
        /**
         * Shell occupation
         */
        occupation?: number;
      }[];
      /**
       * location of the pseudopotential file on filesystem
       */
      path: string;
      /**
       * The names of the simulation engines that can use this pseudopotential, e.g. espresso
       */
      apps: string[];
      /**
       * filename of pseudopotential file on filesystem
       */
      filename?: string;
      /**
       * name of the data category
       */
      name?: "pseudopotential";
    };
    /**
     * TODO: remove in the future
     */
    source?: {
      info?: {};
      type?: string;
    };
  }[];
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/methods_directory/physical/pw.json */

/**
 * Approximating the electronic wave function with a plane wave basis
 */
export interface UnitMethodPlaneWave {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pw";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "wf";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/methods_directory/physical/smearing.json */

/**
 * Approximating Heaviside step function with smooth function
 */
export interface UnitMethodSmearing {
  /**
   * Approximating Heaviside step function with smooth function
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "smearing";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      ("gaussian" | "marzari-vanderbilt" | "methfessel-paxton" | "fermi-dirac");
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "wf";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/methods_directory/physical/tetrahedron.json */

export interface UnitMethodTetrahedron {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "tetrahedron";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      ("linear" | "optimized" | "bloechl");
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "wf";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Instructive parameters defining the method
   */
  parameters?: {};
  /**
   * Object showing the actual possible precision based on theory and implementation
   */
  precision?: {};
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/model/categorized_model.json */

export interface CategorizedModel {
  method: {
    units: {
      /**
       * Used to categorize entities such as models and methods
       */
      categories?: {
        /**
         * contains either object with slugified entry or slug only as a string
         */
        tier1?:
          | {
              /**
               * descriptive human-readable name of entry
               */
              name: string;
              /**
               * machine-readable identifier
               */
              slug: string;
            }
          | string;
        /**
         * contains either object with slugified entry or slug only as a string
         */
        tier2?:
          | {
              /**
               * descriptive human-readable name of entry
               */
              name: string;
              /**
               * machine-readable identifier
               */
              slug: string;
            }
          | string;
        /**
         * contains either object with slugified entry or slug only as a string
         */
        tier3?:
          | {
              /**
               * descriptive human-readable name of entry
               */
              name: string;
              /**
               * machine-readable identifier
               */
              slug: string;
            }
          | string;
        /**
         * contains either object with slugified entry or slug only as a string
         */
        type?:
          | {
              /**
               * descriptive human-readable name of entry
               */
              name: string;
              /**
               * machine-readable identifier
               */
              slug: string;
            }
          | string;
        /**
         * contains either object with slugified entry or slug only as a string
         */
        subtype?:
          | {
              /**
               * descriptive human-readable name of entry
               */
              name: string;
              /**
               * machine-readable identifier
               */
              slug: string;
            }
          | string;
      };
      /**
       * Instructive parameters defining the method
       */
      parameters?: {};
      /**
       * Object showing the actual possible precision based on theory and implementation
       */
      precision?: {};
      /**
       * entity name
       */
      name?: string;
      /**
       * TODO: Use regex once schema draft version has been updated
       */
      path?: string;
      /**
       * entity tags
       */
      tags?: string[];
    }[];
    /**
     * entity name
     */
    name?: string;
    /**
     * TODO: Use regex once schema draft version has been updated
     */
    path?: string;
    /**
     * entity tags
     */
    tags?: string[];
  };
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters: {};
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/model/mixins/dft/double_hybrid_functional.json */

export interface DoubleHybridFunctionalMixin {
  functional?: "b2plyp";
}
 
/** Schema dist/js/schema/model/mixins/dft/enum_options.json */

export interface ModelMixinsDftEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/mixins/dft/gga_functional.json */

export interface GGAFunctionalMixin {
  functional?: "pbe" | "pbesol";
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/mixins/dft/hybrid_functional.json */

export interface HybridFunctionalMixin {
  functional?: "hse06" | "b3lyp";
}
 
/** Schema dist/js/schema/model/mixins/dft/lda_functional.json */

export interface LDAFunctionalMixin {
  functional?: "pz";
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/mixins/dft/mgga_functional.json */

export interface MetaGGAFunctionalMixin {
  functional?: "scan";
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/mixins/dispersion_correction.json */

export interface DispersionCorrectionMixin {
  dispersionCorrection?: "dft-d2" | "dft-d3" | "xdm" | "ts";
}
 
/** Schema dist/js/schema/model/mixins/enum_options.json */

export interface ModelMixinsEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/mixins/hubbard.json */

export interface HubbardModelMixin {
  hubbardType?: "u";
}
 
/** Schema dist/js/schema/model/mixins/spin_orbit_coupling.json */

export interface SpinOrbitCouplingMixin {
  spinOrbitCoupling?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/mixins/spin_polarization.json */

export interface SpinPolarizationMixin {
  spinPolarization?: "collinear" | "non-collinear";
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/model_parameters.json */

export type ModelParameters = ModelParameters1 & ModelParameters2;
export type ModelParameters2 =
  | {
      functional?: "pz";
      [k: string]: unknown;
    }
  | {
      functional?: "pbe" | "pbesol";
      [k: string]: unknown;
    }
  | {
      functional?: "scan";
      [k: string]: unknown;
    }
  | {
      functional?: "hse06" | "b3lyp";
    }
  | {
      functional?: "b2plyp";
    };

export interface ModelParameters1 {
  hubbardType?: "u";
  spinPolarization?: "collinear" | "non-collinear";
  spinOrbitCoupling?: boolean;
  dispersionCorrection?: "dft-d2" | "dft-d3" | "xdm" | "ts";
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/model/model_without_method.json */

export interface ModelWithoutMethodSchemaBase {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters: {};
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/model.json */

export interface BaseModel {
  /**
   * general type of the model, eg. `dft`
   */
  type: string;
  /**
   * general subtype of the model, eg. `lda`
   */
  subtype: string;
  method: {
    /**
     * general type of this method, eg. `pseudopotential`
     */
    type: string;
    /**
     * general subtype of this method, eg. `ultra-soft`
     */
    subtype: string;
    /**
     * Object showing the actual possible precision based on theory and implementation
     */
    precision?: {};
    /**
     * additional data specific to method, eg. array of pseudopotentials
     */
    data?: {};
  };
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/enum_options.json */

export interface ModelsCategoryEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/pb/enum_options.json */

export interface ModelsCategoryPbEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/pb/qm/abin/enum_options.json */

export interface ModelsCategoryPbQmAbinEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/pb/qm/abin/gw.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface GWCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "gw";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    ("g0w0" | "evgw0" | "evgw");
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "abin";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
}
 
/** Schema dist/js/schema/models_category/pb/qm/abin.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface AbInitioCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "abin";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/enum_options.json */

export interface ModelsCategoryPbQmDftEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/ksdft/double_hybrid.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DFTDoubleHybridFunctionalCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "double-hybrid";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ksdft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/ksdft/enum_options.json */

export interface ModelsCategoryPbQmDftKsdftEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/ksdft/gga.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DFTGGAFunctionalCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "gga";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ksdft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/ksdft/hybrid.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DFTHybridFunctionalCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "hybrid";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ksdft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/ksdft/lda.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DFTLDAFunctionalCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "lda";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ksdft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/ksdft/mgga.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DFTMetaGGAFunctionalCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "mgga";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ksdft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft/ksdft.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface KohnShamDFTCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ksdft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/pb/qm/dft.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DensityFunctionalTheoryCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "dft";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/pb/qm/enum_options.json */

export interface ModelsCategoryPbQmEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/pb/qm/semp.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface SemiEmpiricalCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "semp";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/pb/qm.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface QuantumMechanicalCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "qm";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/pb.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface PhysicsBasedModelCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "pb";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/st/det/enum_options.json */

export interface ModelsCategoryStDetEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/st/det/ml/enum_options.json */

export interface ModelsCategoryStDetMlEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/st/det/ml/re.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface RegressionModelCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "re";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ml";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "det";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "st";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/st/det/ml.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface MachineLearningModelCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "ml";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "det";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "st";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/st/det.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface DeterministicModelCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "det";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "st";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_category/st/enum_options.json */

export interface ModelsCategoryStEnumOptions {
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_category/st.json */

/**
 * Used to categorize entities such as models and methods
 */
export interface StatisticalModelCategorySchema {
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier1?: (
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string
  ) &
    "st";
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier2?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  tier3?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  type?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
  /**
   * contains either object with slugified entry or slug only as a string
   */
  subtype?:
    | {
        /**
         * descriptive human-readable name of entry
         */
        name: string;
        /**
         * machine-readable identifier
         */
        slug: string;
      }
    | string;
}
 
/** Schema dist/js/schema/models_directory/double_hybrid.json */

export interface ModelDoubleHybridFunctional {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "double-hybrid";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ksdft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "dft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pb";
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters:
    | {
        spinOrbitCoupling?: boolean;
        [k: string]: unknown;
      }
    | {
        dispersionCorrection?: "dft-d2" | "dft-d3" | "xdm" | "ts";
      }
    | {
        spinPolarization?: "collinear" | "non-collinear";
        [k: string]: unknown;
      };
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/models_directory/gga.json */

export interface ModelGeneralizedGradientApproximation {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "gga";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ksdft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "dft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pb";
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters:
    | {
        spinOrbitCoupling?: boolean;
        [k: string]: unknown;
      }
    | {
        dispersionCorrection?: "dft-d2" | "dft-d3" | "xdm" | "ts";
      }
    | {
        spinPolarization?: "collinear" | "non-collinear";
        [k: string]: unknown;
      }
    | {
        hubbardType?: "u";
      };
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/models_directory/gw.json */

export interface ModelGwApproximation {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "gw";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      ("g0w0" | "evgw0" | "evgw");
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "abin";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pb";
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters: (
    | {
        spinPolarization?: "collinear" | "non-collinear";
        [k: string]: unknown;
      }
    | {
        spinOrbitCoupling?: boolean;
        [k: string]: unknown;
      }
  ) &
    (
      | {
          functional?: "pz";
          [k: string]: unknown;
        }
      | {
          functional?: "pbe" | "pbesol";
          [k: string]: unknown;
        }
      | {
          functional?: "scan";
          [k: string]: unknown;
        }
    );
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/models_directory/hybrid.json */

export interface ModelHybridFunctional {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "hybrid";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ksdft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "dft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pb";
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters:
    | {
        spinOrbitCoupling?: boolean;
        [k: string]: unknown;
      }
    | {
        dispersionCorrection?: "dft-d2" | "dft-d3" | "xdm" | "ts";
      }
    | {
        spinPolarization?: "collinear" | "non-collinear";
        [k: string]: unknown;
      }
    | {
        hubbardType?: "u";
      };
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/models_directory/lda.json */

export interface ModelLocalDensityApproximation {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "lda";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ksdft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "dft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pb";
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters:
    | {
        spinOrbitCoupling?: boolean;
        [k: string]: unknown;
      }
    | {
        dispersionCorrection?: "dft-d2" | "dft-d3" | "xdm" | "ts";
      }
    | {
        spinPolarization?: "collinear" | "non-collinear";
        [k: string]: unknown;
      }
    | {
        hubbardType?: "u";
      };
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/models_directory/legacy/dft.json */

export type LegacyModelDensityFunctionalTheory = LegacyModelDensityFunctionalTheory1 &
  LegacyModelDensityFunctionalTheory2;
export type LegacyModelDensityFunctionalTheory2 =
  | {
      subtype?: "lda";
      functional?: "pz" | "pw" | "vwn" | "other";
    }
  | {
      subtype?: "gga";
      functional?: "pbe" | "pbesol" | "pw91" | "other";
    }
  | {
      subtype?: "hybrid";
      functional?: "b3lyp" | "hse06";
    };

export interface LegacyModelDensityFunctionalTheory1 {
  /**
   * general type of the model, eg. `dft`
   */
  type: "dft";
  /**
   * general subtype of the model, eg. `lda`
   */
  subtype: string;
  method: {
    /**
     * general type of this method, eg. `pseudopotential`
     */
    type: string;
    /**
     * general subtype of this method, eg. `ultra-soft`
     */
    subtype: string;
    /**
     * Object showing the actual possible precision based on theory and implementation
     */
    precision?: {};
    /**
     * additional data specific to method, eg. array of pseudopotentials
     */
    data?: {};
  };
  [k: string]: unknown;
}
/**
 * This interface was referenced by `LegacyModelDensityFunctionalTheory1`'s JSON-Schema
 * via the `definition` "lda".
 */
export interface Lda {
  subtype?: "lda";
  functional?: "pz" | "pw" | "vwn" | "other";
}
/**
 * This interface was referenced by `LegacyModelDensityFunctionalTheory1`'s JSON-Schema
 * via the `definition` "gga".
 */
export interface Gga {
  subtype?: "gga";
  functional?: "pbe" | "pbesol" | "pw91" | "other";
}
/**
 * This interface was referenced by `LegacyModelDensityFunctionalTheory1`'s JSON-Schema
 * via the `definition` "hybrid".
 */
export interface Hybrid {
  subtype?: "hybrid";
  functional?: "b3lyp" | "hse06";
}
 
/** Schema dist/js/schema/models_directory/legacy/ml.json */

export interface LegacyModelRegression {
  /**
   * general type of the model, eg. `dft`
   */
  type: "ml";
  /**
   * general subtype of the model, eg. `lda`
   */
  subtype: "re";
  method: {
    /**
     * general type of this method, eg. `pseudopotential`
     */
    type: string;
    /**
     * general subtype of this method, eg. `ultra-soft`
     */
    subtype: string;
    /**
     * Object showing the actual possible precision based on theory and implementation
     */
    precision?: {};
    /**
     * additional data specific to method, eg. array of pseudopotentials
     */
    data?: {};
  };
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_directory/legacy/unknown.json */

export interface LegacyModelUnknown {
  /**
   * general type of the model, eg. `dft`
   */
  type: "unknown";
  /**
   * general subtype of the model, eg. `lda`
   */
  subtype: "unknown";
  method: {
    /**
     * general type of this method, eg. `pseudopotential`
     */
    type: string;
    /**
     * general subtype of this method, eg. `ultra-soft`
     */
    subtype: string;
    /**
     * Object showing the actual possible precision based on theory and implementation
     */
    precision?: {};
    /**
     * additional data specific to method, eg. array of pseudopotentials
     */
    data?: {};
  };
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/models_directory/mgga.json */

export interface ModelMetaGeneralizedGradientApproximation {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "mgga";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ksdft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "dft";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "qm";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "pb";
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters:
    | {
        spinOrbitCoupling?: boolean;
        [k: string]: unknown;
      }
    | {
        dispersionCorrection?: "dft-d2" | "dft-d3" | "xdm" | "ts";
      }
    | {
        spinPolarization?: "collinear" | "non-collinear";
        [k: string]: unknown;
      }
    | {
        hubbardType?: "u";
      };
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/models_directory/re.json */

/**
 * machine learning model type/subtype schema
 */
export interface ModelRegression {
  /**
   * Used to categorize entities such as models and methods
   */
  categories: {
    /**
     * contains either object with slugified entry or slug only as a string
     */
    type?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "re";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier3?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "ml";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier2?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "det";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    tier1?: (
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string
    ) &
      "st";
    /**
     * contains either object with slugified entry or slug only as a string
     */
    subtype?:
      | {
          /**
           * descriptive human-readable name of entry
           */
          name: string;
          /**
           * machine-readable identifier
           */
          slug: string;
        }
      | string;
  };
  /**
   * Model parameters defined in-place or via model mixins
   */
  parameters: {};
  reference?: {
    type?: "literature";
    /**
     * Digital Object Identifier of the reference.
     */
    doi?: string;
    /**
     * International Standard Book Number of the reference.
     */
    isbn?: string;
    /**
     * International Standard Serial Number of the reference.
     */
    issn?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    /**
     * Publisher of the work.
     */
    publisher?: string;
    /**
     * Journal in which the work appeared.
     */
    journal?: string;
    /**
     * Volume of the series in which the work appeared.
     */
    volume?: string;
    /**
     * Year in which the reference was published.
     */
    year?: string;
    /**
     * Issue of the collection in which the work appeared.
     */
    issue?: string;
    /**
     * Start and end pages of the work.
     */
    pages?: {
      start: string;
      end?: string;
    };
    /**
     * List of authors of the work.
     */
    authors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * List of editors of the work.
     */
    editors?: {
      first: string;
      middle?: string;
      last: string;
      affiliation?: string;
    }[];
    /**
     * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
     */
    reference?: {}[];
  };
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/project.json */

export interface ProjectSchema {
  /**
   * project GID
   */
  gid?: number;
  /**
   * charge rates info for project
   */
  clusterBasedChargeRates?: {
    rate?: number;
    timestamp?: number;
    hostname?: string;
  }[];
  isExternal?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  metadata?: {};
}
 
/** Schema dist/js/schema/properties_directory/derived_properties.json */

export type DerivedPropertiesSchema = (
  | {
      name?: "volume";
      units?: "angstrom^3";
      value: number;
    }
  | {
      name?: "density";
      units?: "g/cm^3";
      value: number;
    }
  | {
      /**
       * point group symbol in Schoenflies notation
       */
      pointGroupSymbol?: string;
      /**
       * space group symbol in Hermann–Mauguin notation
       */
      spaceGroupSymbol?: string;
      /**
       * tolerance used for symmetry calculation
       */
      tolerance?: {
        units?: "angstrom";
        value: number;
      };
      name?: "symmetry";
    }
  | {
      name?: "elemental_ratio";
      value: number;
      /**
       * the element this ratio is for
       */
      element?: string;
    }
  | {
      name?: "p-norm";
      /**
       * degree of the dimensionality of the norm
       */
      degree?: number;
      value: number;
    }
  | {
      name?: "inchi";
      value: string;
    }
  | {
      name?: "inchi_key";
      value: string;
    }
)[];
 
/** Schema dist/js/schema/properties_directory/electronic_configuration.json */

export interface ElectronicConfigurationSchema {
  /**
   * total charge of the molecular system
   */
  charge?: number;
  /**
   * calculated as 2S+1, with S is the total spin angular momentum
   */
  multiplicity?: number;
}
 
/** Schema dist/js/schema/properties_directory/elemental/atomic_radius.json */

/**
 * atomic radius
 */
export interface AtomicRadius {
  name?: "atomic_radius";
  units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/elemental/electronegativity.json */

/**
 * electronegativity for the element (Pauling scale)
 */
export interface Electronegativity {
  name?: "electronegativity";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/elemental/ionization_potential.json */

/**
 * ionization potential for the element
 */
export interface IonizationPotential {
  name?: "ionization_potential";
  units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/average_potential_profile.json */

export interface AveragePotentialProfileSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "z coordinate";
    /**
     * units for an axis
     */
    units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "energy";
    /**
     * units for an axis
     */
    units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
  };
  name?: "average_potential_profile";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/band_gaps.json */

/**
 * contains band gap values
 */
export interface BandGapsSchema {
  name: "band_gaps";
  values?: {
    /**
     * @minItems 3
     * @maxItems 3
     */
    kpointConduction?: [number, number, number];
    /**
     * @minItems 3
     * @maxItems 3
     */
    kpointValence?: [number, number, number];
    /**
     * eigenvalue at k-point in conduction band
     */
    eigenvalueConduction?: number;
    /**
     * eigenvalue at k-point in valence band
     */
    eigenvalueValence?: number;
    spin?: number;
    type: "direct" | "indirect";
    units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
    value: number;
  }[];
  eigenvalues?: {
    /**
     * @minItems 3
     * @maxItems 3
     */
    kpoint?: [number, number, number];
    weight?: number;
    eigenvalues?: {
      spin?: number;
      energies?: unknown[];
      occupations?: unknown[];
    }[];
  }[];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/band_structure.json */

export interface BandStructureSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "kpoints";
    /**
     * units for an axis
     */
    units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "energy";
    /**
     * units for an axis
     */
    units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
  };
  name?: "band_structure";
  /**
   * spin of each band
   */
  spin?: (0.5 | -0.5)[];
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/charge_density_profile.json */

export interface ChargeDensityProfileSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "z coordinate";
    /**
     * units for an axis
     */
    units?: string;
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "charge density";
    /**
     * units for an axis
     */
    units?: "e/A";
  };
  name?: "charge_density_profile";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/density_of_states.json */

export interface DensityOfStatesSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "energy";
    /**
     * units for an axis
     */
    units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "density of states";
    /**
     * units for an axis
     */
    units?: "states/unitcell";
  };
  name?: "density_of_states";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [
    {
      /**
       * chemical element
       */
      element?: string;
      /**
       * index inside sub-array of atoms of the same element type
       */
      index?: number;
      /**
       * electronic character and shell of PDOS, such as `1s` or `s`, or `total`
       */
      electronicState?: string;
      /**
       * spin of the electronic state
       */
      spin?: 0.5 | -0.5;
    },
    ...{
      /**
       * chemical element
       */
      element?: string;
      /**
       * index inside sub-array of atoms of the same element type
       */
      index?: number;
      /**
       * electronic character and shell of PDOS, such as `1s` or `s`, or `total`
       */
      electronicState?: string;
      /**
       * spin of the electronic state
       */
      spin?: 0.5 | -0.5;
    }[]
  ];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/dielectric_tensor.json */

/**
 * The real and imaginary parts of the diagonal elements of the dieletric tensor
 */
export interface DielectricTensorProperty {
  name: "dielectric_tensor";
  values?: {
    /**
     * Real or imaginary part of the dielectric tensor component
     */
    part: "real" | "imaginary";
    spin?: number;
    /**
     * Frequencies
     */
    frequencies: number[];
    /**
     * Matrix with 3 columns, e.g. x, y, z
     */
    components: [number, number, number][];
  }[];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/file_content.json */

export interface FileContent {
  name: "file_content";
  /**
   * What kind of file this is, e.g. image / text
   */
  filetype?: "image" | "text" | "csv";
  objectData: {
    /**
     * Object storage container for the file
     */
    CONTAINER?: string;
    /**
     * Name of the file inside the object storage bucket
     */
    NAME?: string;
    /**
     * Object storage provider
     */
    PROVIDER?: string;
    /**
     * Region for the object container specified in Container
     */
    REGION?: string;
    /**
     * Size of the file in bytes
     */
    SIZE?: number;
    /**
     * Unix timestamp showing when the file was last modified
     */
    TIMESTAMP?: string;
  };
  /**
   * Relative path to the directory that contains the file.
   */
  pathname?: string;
  /**
   * Basename of the file
   */
  basename?: string;
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/hubbard_u.json */

/**
 * Hubbard U values in eV corresponding to atomic species, orbital and site number.
 */
export interface HubbardUParameters {
  name: "hubbard_u";
  units?: "eV";
  values?: {
    /**
     * Site number or index in the lattice
     */
    id: number;
    /**
     * Example: Co1, Mn
     */
    atomicSpecies: string;
    orbitalName: string;
    /**
     * Value related to a specific property, e.g., Hubbard U, V etc.
     */
    value: number;
  }[];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/hubbard_v.json */

/**
 * Hubbard V values corresponding to atomic pairs
 */
export interface HubbardVParameters {
  name: "hubbard_v" | "hubbard_v_nn";
  units?: "eV";
  values?: {
    /**
     * Site number or index in the lattice
     */
    id: number;
    /**
     * Site number or index in the lattice of second site
     */
    id2: number;
    /**
     * Example: Co1, Mn
     */
    atomicSpecies: string;
    /**
     * Example: Co2, O
     */
    atomicSpecies2: string;
    orbitalName?: string;
    orbitalName2?: string;
    /**
     * Distance between two sites in Bohr.
     */
    distance?: number;
    /**
     * Value related to a specific property, e.g., Hubbard U, V etc.
     */
    value: number;
  }[];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/hubbard_v_nn.json */

/**
 * Hubbard V value in eV for nearest neighbors used in hp.x output parsing
 */
export interface HubbardVParametersForNearestNeighbors {
  name: "hubbard_v" | "hubbard_v_nn";
  units?: "eV";
  values?: {
    /**
     * Site number or index in the lattice
     */
    id: number;
    /**
     * Site number or index in the lattice of second site
     */
    id2: number;
    /**
     * Example: Co1, Mn
     */
    atomicSpecies: string;
    /**
     * Example: Co2, O
     */
    atomicSpecies2: string;
    orbitalName?: string;
    orbitalName2?: string;
    /**
     * Distance between two sites in Bohr.
     */
    distance?: number;
    /**
     * Value related to a specific property, e.g., Hubbard U, V etc.
     */
    value: number;
  }[];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/phonon_dispersions.json */

export interface PhononBandStructureSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "qpoints";
    /**
     * units for an axis
     */
    units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "frequency";
    /**
     * units for an axis
     */
    units?: "cm-1" | "THz" | "meV";
  };
  name?: "phonon_dispersions";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/phonon_dos.json */

export interface PhononDensityOfStatesSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "frequency";
    /**
     * units for an axis
     */
    units?: "cm-1" | "THz" | "meV";
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "Phonon DOS";
    /**
     * units for an axis
     */
    units?: "states/cm-1" | "states/THz" | "states/meV";
  };
  name?: "phonon_dos";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/potential_profile.json */

export interface PotentialProfileSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "z coordinate";
    /**
     * units for an axis
     */
    units?: string;
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "energy";
    /**
     * units for an axis
     */
    units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
  };
  name?: "potential_profile";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/reaction_energy_profile.json */

export interface ReactionEnergyProfileSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "reaction coordinate";
    /**
     * units for an axis
     */
    units?: string;
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "energy";
    /**
     * units for an axis
     */
    units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
  };
  name?: "reaction_energy_profile";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/stress_tensor.json */

export interface StressTensorSchema {
  /**
   * @minItems 3
   * @maxItems 3
   */
  value?: [[number, number, number], [number, number, number], [number, number, number]];
  name?: "stress_tensor";
  units?: "kbar" | "pa";
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/total_energy_contributions.json */

export interface TotalEnergyContributionsSchema {
  /**
   * product of temperature and configurational entropy
   */
  temperatureEntropy?: {
    name?: "temperature_entropy";
    value: number;
  };
  /**
   * non self-consitent energy based on an input charge density
   */
  harrisFoulkes?: {
    name?: "harris_foulkes";
    value: number;
  };
  /**
   * kinetic + pseudopotential energy
   */
  oneElectron?: {
    name?: "one_electron";
    value: number;
  };
  /**
   * energy due to coulomb potential
   */
  hartree?: {
    name?: "hartree";
    value: number;
  };
  /**
   * exchange energy
   */
  exchange?: {
    name?: "exchange";
    value: number;
  };
  /**
   * exchange and correlation energy per particle
   */
  exchangeCorrelation?: {
    name?: "exchange_correlation";
    value: number;
  };
  /**
   * summation of interaction energies at long length scales due to coloumbic interactions
   */
  ewald?: {
    name?: "ewald";
    value: number;
  };
  /**
   * divergent electrostatic ion interaction in compensating electron gas
   */
  alphaZ?: {
    name?: "alphaZ";
    value: number;
  };
  /**
   * kinetic energy of wavefunctions in the atomic limit
   */
  atomicEnergy?: {
    name?: "atomic_energy";
    value: number;
  };
  /**
   * sum of one electron energies of kinetic, electrostatic, and exchange correlation
   */
  eigenvalues?: {
    name?: "eigenvalues";
    value: number;
  };
  /**
   * double counting correction 2
   */
  PAWDoubleCounting2?: {
    name?: "PAW_double-counting_correction_2";
    value: number;
  };
  /**
   * double counting correction 3
   */
  PAWDoubleCounting3?: {
    name?: "PAW_double-counting_correction_3";
    value: number;
  };
  /**
   * hartree-fock contribution
   */
  hartreeFock?: {
    name?: "hartree_fock";
    value: number;
  };
  name?: "total_energy_contributions";
  units?: "kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom";
}
 
/** Schema dist/js/schema/properties_directory/non_scalar/vibrational_spectrum.json */

export interface VibrationalSpectrumSchema {
  xAxis: {
    /**
     * label of an axis object
     */
    label: "frequency" | "wavenumber";
    /**
     * units for an axis
     */
    units?: "cm-1" | "THz" | "meV";
  };
  yAxis: {
    /**
     * label of an axis object
     */
    label: "Intensity" | "Absorbance" | "Absorption coefficient";
    /**
     * units for an axis
     */
    units?: "(debye/angstrom)^2" | "km/mol" | "m/mol" | "a.u.";
  };
  name?: "vibrational_spectrum";
  /**
   * Legend of y Axis data series
   *
   * @minItems 1
   */
  legend?: [unknown, ...unknown[]];
  /**
   * array containing values of x Axis
   */
  xDataArray: unknown[];
  yDataSeries: [number | string, ...(number | string)[]][];
}
 
/** Schema dist/js/schema/properties_directory/scalar/electron_affinity.json */

export interface ElectronAffinitySchema {
  name: "electron_affinity";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/fermi_energy.json */

export interface FermiEnergySchema {
  name: "fermi_energy";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/formation_energy.json */

export interface FormationEnergySchema {
  name: "formation_energy";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/ionization_potential.json */

export interface IonizationPotentialSchema {
  name: "ionization_potential";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/pressure.json */

/**
 * average pressure in unit cell
 */
export interface Pressure {
  name?: "pressure";
  units?: "kbar" | "pa";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/reaction_energy_barrier.json */

export interface ReactionEnergyBarrierSchema {
  name: "reaction_energy_barrier";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/surface_energy.json */

export interface SurfaceEnergySchema {
  name: "surface_energy";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/total_energy.json */

export interface TotalEnergySchema {
  name: "total_energy";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/total_force.json */

export interface TotalForcesSchema {
  name?: "total_force";
  units?: "eV/bohr" | "eV/angstrom" | "rydberg/a.u." | "newton" | "kg*m/s^2" | "eV/a.u.";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/valence_band_offset.json */

export interface ValenceBandOffsetSchema {
  name: "valence_band_offset";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/scalar/zero_point_energy.json */

export interface ZeroPointEnergySchema {
  name: "zero_point_energy";
  units: ("kJ/mol" | "eV" | "J/mol" | "hartree" | "cm-1" | "rydberg" | "eV/atom") | "eV/A^2";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/structural/atomic_forces.json */

/**
 * coordinates of atoms by ids, vector, unitless
 */
export interface AtomicForces {
  name?: "atomic_forces";
  /**
   * array of objects containing integer id each
   */
  values?: {
    value?: [number, number, number] | [boolean, boolean, boolean];
    /**
     * integer id of this entry
     */
    id?: number;
  }[];
  units?: "eV/bohr" | "eV/angstrom" | "rydberg/a.u." | "newton" | "kg*m/s^2" | "eV/a.u.";
}
 
/** Schema dist/js/schema/properties_directory/structural/basis/atomic_constraints.json */

/**
 * atomic constraints schema
 */
export interface AtomicConstraints {
  name?: "atomic_constraints";
  /**
   * array of objects containing integer id each
   */
  values?: {
    value?: [number, number, number] | [boolean, boolean, boolean];
    /**
     * integer id of this entry
     */
    id?: number;
  }[];
}
 
/** Schema dist/js/schema/properties_directory/structural/basis/atomic_coordinate.json */

/**
 * coordinates of atoms by ids, vector, unitless
 */
export interface AtomicCoordinate {
  id?: number;
  value?: [number, number, number] | [boolean, boolean, boolean];
}
 
/** Schema dist/js/schema/properties_directory/structural/basis/atomic_coordinates.json */

/**
 * coordinates of atoms by ids, vector, unitless
 */
export interface AtomicCoordinates {
  name?: "atomic_coordinates";
  values?: {
    id?: number;
    value?: [number, number, number] | [boolean, boolean, boolean];
  }[];
  units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
}
 
/** Schema dist/js/schema/properties_directory/structural/basis/atomic_element.json */

/**
 * elements of atoms by ids, string, unitless
 */
export interface AtomicElements {
  id: number;
  value: string;
  /**
   * Occurrence is for fractional occupations
   */
  occurrence?: number;
  oxidationState?: number;
}
 
/** Schema dist/js/schema/properties_directory/structural/basis/bonds.json */

export type BondsSchema = {
  /**
   * indices of the two connected atoms
   *
   * @minItems 2
   * @maxItems 2
   */
  atomPair?: [
    {
      /**
       * integer id of this entry
       */
      id?: number;
    },
    {
      /**
       * integer id of this entry
       */
      id?: number;
    }
  ];
  bondType?: "single" | "double" | "triple" | "quadruple" | "aromatic" | "tautomeric" | "dative" | "other";
}[];
 
/** Schema dist/js/schema/properties_directory/structural/basis.json */

export interface BasisSchema {
  elements: {
    id: number;
    value: string;
    /**
     * Occurrence is for fractional occupations
     */
    occurrence?: number;
    oxidationState?: number;
  }[];
  /**
   * Optional numeric label (e.g., 1, 2, as in Fe1, Fe2) to distinguish same atomic species to attach different spin magnetic moment.
   */
  labels?: {
    id?: number;
    value?: number;
  }[];
  coordinates: {
    id?: number;
    value?: [number, number, number] | [boolean, boolean, boolean];
  }[];
  name?: string;
  units?: string;
  bonds?: {
    /**
     * indices of the two connected atoms
     *
     * @minItems 2
     * @maxItems 2
     */
    atomPair?: [
      {
        /**
         * integer id of this entry
         */
        id?: number;
      },
      {
        /**
         * integer id of this entry
         */
        id?: number;
      }
    ];
    bondType?: "single" | "double" | "triple" | "quadruple" | "aromatic" | "tautomeric" | "dative" | "other";
  }[];
}
 
/** Schema dist/js/schema/properties_directory/structural/density.json */

export interface DensitySchema {
  name?: "density";
  units?: "g/cm^3";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/structural/elemental_ratio.json */

/**
 * ration of this element in the compound
 */
export interface ElementalRatio {
  name?: "elemental_ratio";
  value: number;
  /**
   * the element this ratio is for
   */
  element?: string;
}
 
/** Schema dist/js/schema/properties_directory/structural/inchi.json */

export interface InChIRepresentationSchema {
  name?: "inchi";
  value: string;
}
 
/** Schema dist/js/schema/properties_directory/structural/inchi_key.json */

export interface InChIKeyRepresentationSchema {
  name?: "inchi_key";
  value: string;
}
 
/** Schema dist/js/schema/properties_directory/structural/lattice/lattice_bravais.json */

export interface LatticeImplicitSchema {
  type:
    | "CUB"
    | "BCC"
    | "FCC"
    | "TET"
    | "MCL"
    | "ORC"
    | "ORCC"
    | "ORCF"
    | "ORCI"
    | "HEX"
    | "BCT"
    | "TRI"
    | "MCLC"
    | "RHL";
  units?: {
    length?: "angstrom" | "bohr";
    angle?: "degree" | "radian";
  };
  /**
   * length of the first lattice vector
   */
  a: number;
  /**
   * length of the second lattice vector
   */
  b: number;
  /**
   * length of the third lattice vector
   */
  c: number;
  /**
   * angle between first and second lattice vector
   */
  alpha: number;
  /**
   * angle between second and third lattice vector
   */
  beta: number;
  /**
   * angle between first and third lattice vector
   */
  gamma: number;
}
 
/** Schema dist/js/schema/properties_directory/structural/lattice/lattice_vectors.json */

export interface LatticeExplicitUnit {
  /**
   * lattice parameter for fractional coordinates
   */
  alat?: number;
  units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
  /**
   * @minItems 3
   * @maxItems 3
   */
  a: [number, number, number];
  /**
   * @minItems 3
   * @maxItems 3
   */
  b: [number, number, number];
  /**
   * @minItems 3
   * @maxItems 3
   */
  c: [number, number, number];
}
 
/** Schema dist/js/schema/properties_directory/structural/lattice/type_enum.json */

export type LatticeTypeSchema =
  | "CUB"
  | "BCC"
  | "FCC"
  | "TET"
  | "MCL"
  | "ORC"
  | "ORCC"
  | "ORCF"
  | "ORCI"
  | "HEX"
  | "BCT"
  | "TRI"
  | "MCLC"
  | "RHL";
 
/** Schema dist/js/schema/properties_directory/structural/lattice/type_extended_enum.json */

export type LatticeTypeExtendedSchema =
  | "BCC"
  | "BCT-1"
  | "BCT-2"
  | "CUB"
  | "FCC"
  | "HEX"
  | "MCL"
  | "MCLC-1"
  | "MCLC-2"
  | "MCLC-3"
  | "MCLC-4"
  | "MCLC-5"
  | "ORC"
  | "ORCC"
  | "ORCF-1"
  | "ORCF-2"
  | "ORCF-3"
  | "ORCI"
  | "RHL-1"
  | "RHL-2"
  | "TET"
  | "TRI_1a"
  | "TRI_2a"
  | "TRI_1b";
 
/** Schema dist/js/schema/properties_directory/structural/lattice.json */

export interface LatticeSchema {
  name?: "lattice";
  vectors?: {
    /**
     * lattice parameter for fractional coordinates
     */
    alat?: number;
    units?: "km" | "m" | "pm" | "nm" | "angstrom" | "a.u." | "bohr" | "fractional" | "crystal" | "cartesian" | "alat";
    /**
     * @minItems 3
     * @maxItems 3
     */
    a: [number, number, number];
    /**
     * @minItems 3
     * @maxItems 3
     */
    b: [number, number, number];
    /**
     * @minItems 3
     * @maxItems 3
     */
    c: [number, number, number];
  };
  type:
    | "CUB"
    | "BCC"
    | "FCC"
    | "TET"
    | "MCL"
    | "ORC"
    | "ORCC"
    | "ORCF"
    | "ORCI"
    | "HEX"
    | "BCT"
    | "TRI"
    | "MCLC"
    | "RHL";
  units?: {
    length?: "angstrom" | "bohr";
    angle?: "degree" | "radian";
  };
  /**
   * length of the first lattice vector
   */
  a: number;
  /**
   * length of the second lattice vector
   */
  b: number;
  /**
   * length of the third lattice vector
   */
  c: number;
  /**
   * angle between first and second lattice vector
   */
  alpha: number;
  /**
   * angle between second and third lattice vector
   */
  beta: number;
  /**
   * angle between first and third lattice vector
   */
  gamma: number;
}
 
/** Schema dist/js/schema/properties_directory/structural/magnetic_moments.json */

/**
 * magnetization on each ion
 */
export interface MagneticMoments {
  name?: "magnetic_moments";
  /**
   * array of objects containing integer id each
   */
  values?: {
    value?: [number, number, number] | [boolean, boolean, boolean];
    /**
     * integer id of this entry
     */
    id?: number;
  }[];
  units?: "uB";
}
 
/** Schema dist/js/schema/properties_directory/structural/molecular_pattern.json */

export type MolecularPatternSchema = (
  | {
      name?: "functional_group";
      /**
       * array of objects containing integer id each
       */
      atoms?: {
        /**
         * whether atom connects to atoms outside of functional group.
         */
        isConnector?: boolean;
        /**
         * integer id of this entry
         */
        id?: number;
      }[];
      /**
       * SMARTS string for classification of FG; https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification
       */
      SMARTS?: string;
    }
  | {
      name?: "ring";
      /**
       * array of objects containing integer id each
       */
      atoms?: {
        /**
         * whether atom connects to atoms outside of functional group.
         */
        isConnector?: boolean;
        /**
         * integer id of this entry
         */
        id?: number;
      }[];
      isAromatic?: boolean;
    }
  | {
      name?: "special_bond";
      /**
       * array of objects containing integer id each
       */
      atoms?: {
        /**
         * whether atom connects to atoms outside of functional group.
         */
        isConnector?: boolean;
        /**
         * integer id of this entry
         */
        id?: number;
      }[];
    }
)[];
 
/** Schema dist/js/schema/properties_directory/structural/p_norm.json */

/**
 * https://en.wikipedia.org/wiki/Norm_(mathematics)#p-norm
 */
export interface PNorm {
  name?: "p-norm";
  /**
   * degree of the dimensionality of the norm
   */
  degree?: number;
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/structural/patterns/functional_group.json */

export interface FunctionalGroupPatternSchema {
  name?: "functional_group";
  /**
   * array of objects containing integer id each
   */
  atoms?: {
    /**
     * whether atom connects to atoms outside of functional group.
     */
    isConnector?: boolean;
    /**
     * integer id of this entry
     */
    id?: number;
  }[];
  /**
   * SMARTS string for classification of FG; https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification
   */
  SMARTS?: string;
}
 
/** Schema dist/js/schema/properties_directory/structural/patterns/ring.json */

export interface RingPatternSchema {
  name?: "ring";
  /**
   * array of objects containing integer id each
   */
  atoms?: {
    /**
     * whether atom connects to atoms outside of functional group.
     */
    isConnector?: boolean;
    /**
     * integer id of this entry
     */
    id?: number;
  }[];
  isAromatic?: boolean;
}
 
/** Schema dist/js/schema/properties_directory/structural/patterns/special_bond.json */

/**
 * Any bonding interaction that cannot be described by simple 2-atom picture, e.g. 3-center-2-electron bond in diborane
 */
export interface SpecialBondPatternSchema {
  name?: "special_bond";
  /**
   * array of objects containing integer id each
   */
  atoms?: {
    /**
     * whether atom connects to atoms outside of functional group.
     */
    isConnector?: boolean;
    /**
     * integer id of this entry
     */
    id?: number;
  }[];
}
 
/** Schema dist/js/schema/properties_directory/structural/symmetry.json */

export interface SymmetrySchema {
  /**
   * point group symbol in Schoenflies notation
   */
  pointGroupSymbol?: string;
  /**
   * space group symbol in Hermann–Mauguin notation
   */
  spaceGroupSymbol?: string;
  /**
   * tolerance used for symmetry calculation
   */
  tolerance?: {
    units?: "angstrom";
    value: number;
  };
  name?: "symmetry";
}
 
/** Schema dist/js/schema/properties_directory/structural/volume.json */

export interface VolumeSchema {
  name?: "volume";
  units?: "angstrom^3";
  value: number;
}
 
/** Schema dist/js/schema/properties_directory/workflow/convergence/electronic.json */

export interface ElectronicSelfConsistencyConvergenceSchema {
  name?: "convergence_electronic";
  units?: "eV" | "rydberg" | "hartree";
  data: number[][];
}
 
/** Schema dist/js/schema/properties_directory/workflow/convergence/ionic.json */

export interface IonicConvergenceSchema {
  name?: "convergence_ionic";
  /**
   * for ionic convergence tolerance shows force tolerance
   */
  tolerance?: {
    [k: string]: unknown;
  };
  /**
   * units for force tolerance
   */
  units?: "eV";
  /**
   * energetic and structural information
   */
  data: {
    /**
     * converged electronic energy for this structure (last in `electronic`)
     */
    energy?: number;
    /**
     * TODO: structural information at each step to be here
     */
    structure?: {};
    /**
     * data about electronic at this ionic step
     */
    electronic?: {
      /**
       * units for force tolerance
       */
      units?: "eV" | "rydberg" | "hartree";
      data?: number[];
    };
  }[];
}
 
/** Schema dist/js/schema/properties_directory/workflow/convergence/kpoint.json */

export interface ConvergenceSchemaForConvergingAPropertyWrtKpoints {
  /**
   * tolerance for the property under investigation
   */
  tolerance: {
    [k: string]: unknown;
  };
  /**
   * units for the property under investigation
   */
  units: string;
  /**
   * name of the property under investigation
   */
  property?: string;
  /**
   * kpoint grid and property information
   */
  data: {
    /**
     * value of the property at this step
     */
    value: {
      [k: string]: unknown;
    };
    /**
     * information about the kpoint grid
     */
    grid: {};
    /**
     * optional kpoint spacing information
     */
    spacing?: number;
  }[];
}
 
/** Schema dist/js/schema/property/base.json */

export interface SchemaOfBaseMaterialSPreliminaryProperty {
  /**
   * property slug, e.g. total_energy
   */
  slug?: string;
  /**
   * property group, e.g. qe:dft:gga:pbe
   */
  group?: string;
  /**
   * container of the information, specific to each property
   */
  data: {};
  source: {
    /**
     * Type of the material property's source.
     */
    type?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    info?:
      | {
          /**
           * Material's identity. Used for protoProperties.
           */
          materialId?: string;
          /**
           * Job's identity
           */
          jobId?: string;
          /**
           * Id of the unit that extracted the result
           */
          unitId?: string;
        }
      | {
          type?: "experiment";
          /**
           * experiment authors
           */
          authors: {
            first: string;
            middle?: string;
            last: string;
            affiliation?: string;
          }[];
          /**
           * method used in experiment
           */
          method: string;
          conditions: {
            /**
             * condition unit
             */
            units?: string;
            /**
             * array of condition values
             */
            scalar?: {
              value?: string;
            }[];
            /**
             * human-readable name of the condition
             */
            name: string;
          }[];
          location?: {
            /**
             * location latitude
             */
            latitude: number;
            /**
             * location longitude
             */
            longitude: number;
          };
          /**
           * epoch time.
           */
          timestamp: number;
          /**
           * Note about experiment
           */
          note?: string;
          /**
           * references to literature articles
           */
          references?: {
            type?: "literature";
            /**
             * Digital Object Identifier of the reference.
             */
            doi?: string;
            /**
             * International Standard Book Number of the reference.
             */
            isbn?: string;
            /**
             * International Standard Serial Number of the reference.
             */
            issn?: string;
            /**
             * Internet address of the reference.
             */
            url?: string;
            /**
             * Publisher of the work.
             */
            publisher?: string;
            /**
             * Journal in which the work appeared.
             */
            journal?: string;
            /**
             * Volume of the series in which the work appeared.
             */
            volume?: string;
            /**
             * Year in which the reference was published.
             */
            year?: string;
            /**
             * Issue of the collection in which the work appeared.
             */
            issue?: string;
            /**
             * Start and end pages of the work.
             */
            pages?: {
              start: string;
              end?: string;
            };
            /**
             * List of authors of the work.
             */
            authors?: {
              first: string;
              middle?: string;
              last: string;
              affiliation?: string;
            }[];
            /**
             * List of editors of the work.
             */
            editors?: {
              first: string;
              middle?: string;
              last: string;
              affiliation?: string;
            }[];
            /**
             * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
             */
            reference?: {}[];
          }[];
        };
  };
  /**
   * Id of the corresponding item in the entity bank that this property is obtained for
   */
  exabyteId?: string[];
  precision?: {};
  /**
   * total number of properties among which this property is the best.
   */
  count?: number;
  /**
   * property system tags, marks property system characteristics, values refined or best (could be both)
   */
  systemTags?: ("isRefined" | "isBest")[];
  /**
   * entity identity
   */
  _id?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
}
 
/** Schema dist/js/schema/property/meta.json */

export interface SchemaOfMaterialSMetaProperties {
  /**
   * property slug, e.g. total_energy
   */
  slug?: string;
  /**
   * property group, e.g. qe:dft:gga:pbe
   */
  group?: string;
  /**
   * container of the information, specific to each property
   */
  data: {};
  source: {
    /**
     * Type of the material property's source.
     */
    type?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    info?:
      | {
          /**
           * Material's identity. Used for protoProperties.
           */
          materialId?: string;
          /**
           * Job's identity
           */
          jobId?: string;
          /**
           * Id of the unit that extracted the result
           */
          unitId?: string;
        }
      | {
          type?: "experiment";
          /**
           * experiment authors
           */
          authors: {
            first: string;
            middle?: string;
            last: string;
            affiliation?: string;
          }[];
          /**
           * method used in experiment
           */
          method: string;
          conditions: {
            /**
             * condition unit
             */
            units?: string;
            /**
             * array of condition values
             */
            scalar?: {
              value?: string;
            }[];
            /**
             * human-readable name of the condition
             */
            name: string;
          }[];
          location?: {
            /**
             * location latitude
             */
            latitude: number;
            /**
             * location longitude
             */
            longitude: number;
          };
          /**
           * epoch time.
           */
          timestamp: number;
          /**
           * Note about experiment
           */
          note?: string;
          /**
           * references to literature articles
           */
          references?: {
            type?: "literature";
            /**
             * Digital Object Identifier of the reference.
             */
            doi?: string;
            /**
             * International Standard Book Number of the reference.
             */
            isbn?: string;
            /**
             * International Standard Serial Number of the reference.
             */
            issn?: string;
            /**
             * Internet address of the reference.
             */
            url?: string;
            /**
             * Publisher of the work.
             */
            publisher?: string;
            /**
             * Journal in which the work appeared.
             */
            journal?: string;
            /**
             * Volume of the series in which the work appeared.
             */
            volume?: string;
            /**
             * Year in which the reference was published.
             */
            year?: string;
            /**
             * Issue of the collection in which the work appeared.
             */
            issue?: string;
            /**
             * Start and end pages of the work.
             */
            pages?: {
              start: string;
              end?: string;
            };
            /**
             * List of authors of the work.
             */
            authors?: {
              first: string;
              middle?: string;
              last: string;
              affiliation?: string;
            }[];
            /**
             * List of editors of the work.
             */
            editors?: {
              first: string;
              middle?: string;
              last: string;
              affiliation?: string;
            }[];
            /**
             * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
             */
            reference?: {}[];
          }[];
        };
  };
  /**
   * Id of the corresponding item in the entity bank that this property is obtained for
   */
  exabyteId?: string[];
  precision?: {};
  /**
   * total number of properties among which this property is the best.
   */
  count?: number;
  /**
   * property system tags, marks property system characteristics, values refined or best (could be both)
   */
  systemTags?: ("isRefined" | "isBest")[];
  /**
   * entity identity
   */
  _id?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
}
 
/** Schema dist/js/schema/property/raw.json */

export interface SchemaOfMaterialSPreliminaryProperty {
  /**
   * property slug, e.g. total_energy
   */
  slug?: string;
  /**
   * property group, e.g. qe:dft:gga:pbe
   */
  group?: string;
  /**
   * container of the information, specific to each property
   */
  data: {};
  source: {
    /**
     * Type of the material property's source.
     */
    type?: string;
    /**
     * Internet address of the reference.
     */
    url?: string;
    info?:
      | {
          /**
           * Material's identity. Used for protoProperties.
           */
          materialId?: string;
          /**
           * Job's identity
           */
          jobId?: string;
          /**
           * Id of the unit that extracted the result
           */
          unitId?: string;
        }
      | {
          type?: "experiment";
          /**
           * experiment authors
           */
          authors: {
            first: string;
            middle?: string;
            last: string;
            affiliation?: string;
          }[];
          /**
           * method used in experiment
           */
          method: string;
          conditions: {
            /**
             * condition unit
             */
            units?: string;
            /**
             * array of condition values
             */
            scalar?: {
              value?: string;
            }[];
            /**
             * human-readable name of the condition
             */
            name: string;
          }[];
          location?: {
            /**
             * location latitude
             */
            latitude: number;
            /**
             * location longitude
             */
            longitude: number;
          };
          /**
           * epoch time.
           */
          timestamp: number;
          /**
           * Note about experiment
           */
          note?: string;
          /**
           * references to literature articles
           */
          references?: {
            type?: "literature";
            /**
             * Digital Object Identifier of the reference.
             */
            doi?: string;
            /**
             * International Standard Book Number of the reference.
             */
            isbn?: string;
            /**
             * International Standard Serial Number of the reference.
             */
            issn?: string;
            /**
             * Internet address of the reference.
             */
            url?: string;
            /**
             * Publisher of the work.
             */
            publisher?: string;
            /**
             * Journal in which the work appeared.
             */
            journal?: string;
            /**
             * Volume of the series in which the work appeared.
             */
            volume?: string;
            /**
             * Year in which the reference was published.
             */
            year?: string;
            /**
             * Issue of the collection in which the work appeared.
             */
            issue?: string;
            /**
             * Start and end pages of the work.
             */
            pages?: {
              start: string;
              end?: string;
            };
            /**
             * List of authors of the work.
             */
            authors?: {
              first: string;
              middle?: string;
              last: string;
              affiliation?: string;
            }[];
            /**
             * List of editors of the work.
             */
            editors?: {
              first: string;
              middle?: string;
              last: string;
              affiliation?: string;
            }[];
            /**
             * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
             */
            reference?: {}[];
          }[];
        };
  };
  /**
   * Id of the corresponding item in the entity bank that this property is obtained for
   */
  exabyteId?: string[];
  precision?: {};
  /**
   * total number of properties among which this property is the best.
   */
  count?: number;
  /**
   * property system tags, marks property system characteristics, values refined or best (could be both)
   */
  systemTags?: ("isRefined" | "isBest")[];
  /**
   * entity identity
   */
  _id?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
}
 
/** Schema dist/js/schema/property/source.json */

export interface TheSourceOfAPropertyThisCouldBeAnArticleASimulationOnExabyteAnExternalSimulationEtc {
  /**
   * Type of the material property's source.
   */
  type?: string;
  /**
   * Internet address of the reference.
   */
  url?: string;
  info?:
    | {
        /**
         * Material's identity. Used for protoProperties.
         */
        materialId?: string;
        /**
         * Job's identity
         */
        jobId?: string;
        /**
         * Id of the unit that extracted the result
         */
        unitId?: string;
      }
    | {
        type?: "experiment";
        /**
         * experiment authors
         */
        authors: {
          first: string;
          middle?: string;
          last: string;
          affiliation?: string;
        }[];
        /**
         * method used in experiment
         */
        method: string;
        conditions: {
          /**
           * condition unit
           */
          units?: string;
          /**
           * array of condition values
           */
          scalar?: {
            value?: string;
          }[];
          /**
           * human-readable name of the condition
           */
          name: string;
        }[];
        location?: {
          /**
           * location latitude
           */
          latitude: number;
          /**
           * location longitude
           */
          longitude: number;
        };
        /**
         * epoch time.
         */
        timestamp: number;
        /**
         * Note about experiment
         */
        note?: string;
        /**
         * references to literature articles
         */
        references?: {
          type?: "literature";
          /**
           * Digital Object Identifier of the reference.
           */
          doi?: string;
          /**
           * International Standard Book Number of the reference.
           */
          isbn?: string;
          /**
           * International Standard Serial Number of the reference.
           */
          issn?: string;
          /**
           * Internet address of the reference.
           */
          url?: string;
          /**
           * Publisher of the work.
           */
          publisher?: string;
          /**
           * Journal in which the work appeared.
           */
          journal?: string;
          /**
           * Volume of the series in which the work appeared.
           */
          volume?: string;
          /**
           * Year in which the reference was published.
           */
          year?: string;
          /**
           * Issue of the collection in which the work appeared.
           */
          issue?: string;
          /**
           * Start and end pages of the work.
           */
          pages?: {
            start: string;
            end?: string;
          };
          /**
           * List of authors of the work.
           */
          authors?: {
            first: string;
            middle?: string;
            last: string;
            affiliation?: string;
          }[];
          /**
           * List of editors of the work.
           */
          editors?: {
            first: string;
            middle?: string;
            last: string;
            affiliation?: string;
          }[];
          /**
           * References cited by the work. Reference objects can nest as deeply as needed. This is useful, for example, when tracking the history of a value referenced in a scholarly article; the top level reference would contain information about where the data was accessed while the nested reference would contain information about where it was originally published.
           */
          reference?: {}[];
        }[];
      };
}
 
/** Schema dist/js/schema/software/application.json */

export interface ApplicationSchemaBase {
  /**
   * The short name of the application. e.g. qe
   */
  shortName?: string;
  /**
   * Application's short description.
   */
  summary?: string;
  /**
   * Application version. e.g. 5.3.5
   */
  version?: string;
  /**
   * Application build. e.g. VTST
   */
  build?: string;
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * Whether licensing is present
   */
  isLicensed?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software/executable.json */

export interface ExecutableSchema {
  /**
   * The name of the executable. e.g. pw.x
   */
  name: string;
  /**
   * _ids of the application this executable belongs to
   */
  applicationId?: string[];
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
}
 
/** Schema dist/js/schema/software/flavor.json */

export interface FlavorSchema {
  /**
   * _id of the executable this flavor belongs to
   */
  executableId?: string;
  /**
   * name of the executable this flavor belongs to
   */
  executableName?: string;
  /**
   * name of the application this flavor belongs to
   */
  applicationName?: string;
  input?: {
    templateId?: string;
    templateName?: string;
    /**
     * name of the resulting input file, if different than template name
     */
    name?: string;
  }[];
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
}
 
/** Schema dist/js/schema/software/template.json */

export interface TemplateSchema {
  applicationName?: string;
  applicationVersion?: string;
  executableName?: string;
  contextProviders?: {
    /**
     * The name of this item. e.g. scf_accuracy
     */
    name: string;
  }[];
  /**
   * Input file name. e.g. pw_scf.in
   */
  name: string;
  /**
   * Content of the input file. e.g. &CONTROL    calculation='scf' ...
   */
  content: string;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
}
 
/** Schema dist/js/schema/software_directory/ml/exabyteml.json */

export interface ExabyteMachineLearningEngineSchema {
  name?: "exabyteml";
  summary?: "exabyte machine learning engine";
  version?: "0.2.0";
}
 
/** Schema dist/js/schema/software_directory/ml/unit/execution/evaluate/cross_validate.json */

export interface CrossValidationUnitSchema {
  /**
   * TODO: consider keeping executable `evaluate` and flavor `cross-validate` as before
   */
  input: {
    /**
     * number of groups to split the training dataset for cross-validation
     */
    nSplits: number;
  };
  /**
   * type of the unit
   */
  type: "execution";
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  executable?: {
    /**
     * The name of the executable. e.g. pw.x
     */
    name: string;
    /**
     * _ids of the application this executable belongs to
     */
    applicationId?: string[];
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  flavor?: {
    /**
     * _id of the executable this flavor belongs to
     */
    executableId?: string;
    /**
     * name of the executable this flavor belongs to
     */
    executableName?: string;
    /**
     * name of the application this flavor belongs to
     */
    applicationName?: string;
    input?: {
      templateId?: string;
      templateName?: string;
      /**
       * name of the resulting input file, if different than template name
       */
      name?: string;
    }[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/ml/unit/execution/initialize.json */

export interface InitializeUnitSchema {
  /**
   * model init unit (NOTE: info about method, eg. regression/linear is taken from (sub)workflow)
   */
  input: {
    /**
     * target properties to predict (NOTE: must be a subset of targets for which training was done)
     */
    targets: string[];
  };
  /**
   * type of the unit
   */
  type: "execution";
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  executable?: {
    /**
     * The name of the executable. e.g. pw.x
     */
    name: string;
    /**
     * _ids of the application this executable belongs to
     */
    applicationId?: string[];
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  flavor?: {
    /**
     * _id of the executable this flavor belongs to
     */
    executableId?: string;
    /**
     * name of the executable this flavor belongs to
     */
    executableName?: string;
    /**
     * name of the application this flavor belongs to
     */
    applicationName?: string;
    input?: {
      templateId?: string;
      templateName?: string;
      /**
       * name of the resulting input file, if different than template name
       */
      name?: string;
    }[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/ml/unit/execution/score.json */

export interface TrainScoreSchema {
  /**
   * unit input (type to be specified by the application's execution unit)
   */
  input: {
    [k: string]: unknown;
  };
  /**
   * type of the unit
   */
  type: "execution";
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  executable?: {
    /**
     * The name of the executable. e.g. pw.x
     */
    name: string;
    /**
     * _ids of the application this executable belongs to
     */
    applicationId?: string[];
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  flavor?: {
    /**
     * _id of the executable this flavor belongs to
     */
    executableId?: string;
    /**
     * name of the executable this flavor belongs to
     */
    executableName?: string;
    /**
     * name of the application this flavor belongs to
     */
    applicationName?: string;
    input?: {
      templateId?: string;
      templateName?: string;
      /**
       * name of the resulting input file, if different than template name
       */
      name?: string;
    }[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/ml/unit/execution/train.json */

export interface TrainUnitSchema {
  /**
   * model train unit (NOTE: info about method, eg. regression/linear is taken from (sub)workflow)
   */
  input: {
    /**
     * material features used for model fitting
     */
    features: string[];
    /**
     * target properties to train for
     */
    targets: string[];
  };
  /**
   * type of the unit
   */
  type: "execution";
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  executable?: {
    /**
     * The name of the executable. e.g. pw.x
     */
    name: string;
    /**
     * _ids of the application this executable belongs to
     */
    applicationId?: string[];
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  flavor?: {
    /**
     * _id of the executable this flavor belongs to
     */
    executableId?: string;
    /**
     * name of the executable this flavor belongs to
     */
    executableName?: string;
    /**
     * name of the application this flavor belongs to
     */
    applicationName?: string;
    input?: {
      templateId?: string;
      templateName?: string;
      /**
       * name of the resulting input file, if different than template name
       */
      name?: string;
    }[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/ml/unit/execution.json */

export type SoftwareDirectoryMlUnitExecution =
  | {
      /**
       * TODO: consider keeping executable `evaluate` and flavor `cross-validate` as before
       */
      input: {
        /**
         * number of groups to split the training dataset for cross-validation
         */
        nSplits: number;
      };
      /**
       * type of the unit
       */
      type: "execution";
      application: {
        /**
         * The short name of the application. e.g. qe
         */
        shortName?: string;
        /**
         * Application's short description.
         */
        summary?: string;
        /**
         * Application version. e.g. 5.3.5
         */
        version?: string;
        /**
         * Application build. e.g. VTST
         */
        build?: string;
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * Whether licensing is present
         */
        isLicensed?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        [k: string]: unknown;
      };
      executable?: {
        /**
         * The name of the executable. e.g. pw.x
         */
        name: string;
        /**
         * _ids of the application this executable belongs to
         */
        applicationId?: string[];
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      flavor?: {
        /**
         * _id of the executable this flavor belongs to
         */
        executableId?: string;
        /**
         * name of the executable this flavor belongs to
         */
        executableName?: string;
        /**
         * name of the application this flavor belongs to
         */
        applicationName?: string;
        input?: {
          templateId?: string;
          templateName?: string;
          /**
           * name of the resulting input file, if different than template name
           */
          name?: string;
        }[];
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * model train unit (NOTE: info about method, eg. regression/linear is taken from (sub)workflow)
       */
      input: {
        /**
         * material features used for model fitting
         */
        features: string[];
        /**
         * target properties to train for
         */
        targets: string[];
      };
      /**
       * type of the unit
       */
      type: "execution";
      application: {
        /**
         * The short name of the application. e.g. qe
         */
        shortName?: string;
        /**
         * Application's short description.
         */
        summary?: string;
        /**
         * Application version. e.g. 5.3.5
         */
        version?: string;
        /**
         * Application build. e.g. VTST
         */
        build?: string;
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * Whether licensing is present
         */
        isLicensed?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        [k: string]: unknown;
      };
      executable?: {
        /**
         * The name of the executable. e.g. pw.x
         */
        name: string;
        /**
         * _ids of the application this executable belongs to
         */
        applicationId?: string[];
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      flavor?: {
        /**
         * _id of the executable this flavor belongs to
         */
        executableId?: string;
        /**
         * name of the executable this flavor belongs to
         */
        executableName?: string;
        /**
         * name of the application this flavor belongs to
         */
        applicationName?: string;
        input?: {
          templateId?: string;
          templateName?: string;
          /**
           * name of the resulting input file, if different than template name
           */
          name?: string;
        }[];
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * unit input (type to be specified by the application's execution unit)
       */
      input: {
        [k: string]: unknown;
      };
      /**
       * type of the unit
       */
      type: "execution";
      application: {
        /**
         * The short name of the application. e.g. qe
         */
        shortName?: string;
        /**
         * Application's short description.
         */
        summary?: string;
        /**
         * Application version. e.g. 5.3.5
         */
        version?: string;
        /**
         * Application build. e.g. VTST
         */
        build?: string;
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * Whether licensing is present
         */
        isLicensed?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        [k: string]: unknown;
      };
      executable?: {
        /**
         * The name of the executable. e.g. pw.x
         */
        name: string;
        /**
         * _ids of the application this executable belongs to
         */
        applicationId?: string[];
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      flavor?: {
        /**
         * _id of the executable this flavor belongs to
         */
        executableId?: string;
        /**
         * name of the executable this flavor belongs to
         */
        executableName?: string;
        /**
         * name of the application this flavor belongs to
         */
        applicationName?: string;
        input?: {
          templateId?: string;
          templateName?: string;
          /**
           * name of the resulting input file, if different than template name
           */
          name?: string;
        }[];
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * model init unit (NOTE: info about method, eg. regression/linear is taken from (sub)workflow)
       */
      input: {
        /**
         * target properties to predict (NOTE: must be a subset of targets for which training was done)
         */
        targets: string[];
      };
      /**
       * type of the unit
       */
      type: "execution";
      application: {
        /**
         * The short name of the application. e.g. qe
         */
        shortName?: string;
        /**
         * Application's short description.
         */
        summary?: string;
        /**
         * Application version. e.g. 5.3.5
         */
        version?: string;
        /**
         * Application build. e.g. VTST
         */
        build?: string;
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * Whether licensing is present
         */
        isLicensed?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        [k: string]: unknown;
      };
      executable?: {
        /**
         * The name of the executable. e.g. pw.x
         */
        name: string;
        /**
         * _ids of the application this executable belongs to
         */
        applicationId?: string[];
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      flavor?: {
        /**
         * _id of the executable this flavor belongs to
         */
        executableId?: string;
        /**
         * name of the executable this flavor belongs to
         */
        executableName?: string;
        /**
         * name of the application this flavor belongs to
         */
        applicationName?: string;
        input?: {
          templateId?: string;
          templateName?: string;
          /**
           * name of the resulting input file, if different than template name
           */
          name?: string;
        }[];
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    };
 
/** Schema dist/js/schema/software_directory/ml/unit/processing/data_transformation/manipulation.json */

export interface ManipulationUnitSchema {
  /**
   * Contains information about the operation used.
   */
  operation: "data_transformation";
  /**
   * Contains information about the specific type of the operation used.
   */
  operationType: "manipulation";
  /**
   * unit input (type to be specified by the child units)
   */
  inputData: {
    /**
     * whether to clean missing data, eg. NaN
     */
    cleanMissingData: boolean;
    /**
     * whether to remove duplicate rows
     */
    removeDuplicateRows: boolean;
    /**
     * replace None values with a given value
     */
    replaceNoneValuesWith: number;
  };
  /**
   * type of the unit
   */
  type: "processing";
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/ml/unit/processing/data_transformation/scale_and_reduce.json */

export interface ScaleAndReduceUnitSchema {
  /**
   * Contains information about the operation used.
   */
  operation: "data_transformation";
  /**
   * Contains information about the specific type of the operation used.
   */
  operationType: "scale_and_reduce";
  /**
   * unit input (type to be specified by the child units)
   */
  inputData: {
    /**
     * type of scaler to be applied
     */
    scaler: "standard_scaler";
    /**
     * per-feature scaling data
     */
    perFeature?: {
      /**
       * variance in original training data
       */
      variance?: number;
      /**
       * mean value of the original training data
       */
      mean?: number;
      /**
       * scale multiplier for this feature/property
       */
      scale: number;
      /**
       * feature/property name in 'flattened' format
       */
      name: string;
    }[];
  };
  /**
   * type of the unit
   */
  type: "processing";
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/ml/unit/processing/data_transformation.json */

export type SoftwareDirectoryMlUnitProcessingDataTransformation = {
  /**
   * Contains information about the operation used.
   */
  operation: "data_transformation";
  /**
   * Contains information about the specific type of the operation used.
   */
  operationType: "scale_and_reduce";
  /**
   * unit input (type to be specified by the child units)
   */
  inputData: {
    /**
     * type of scaler to be applied
     */
    scaler: "standard_scaler";
    /**
     * per-feature scaling data
     */
    perFeature?: {
      /**
       * variance in original training data
       */
      variance?: number;
      /**
       * mean value of the original training data
       */
      mean?: number;
      /**
       * scale multiplier for this feature/property
       */
      scale: number;
      /**
       * feature/property name in 'flattened' format
       */
      name: string;
    }[];
  };
  /**
   * type of the unit
   */
  type: "processing";
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
};
 
/** Schema dist/js/schema/software_directory/ml/unit/processing/feature_selection/filter_based.json */

export interface FilterBasedFeatureSelectionUnitSchema {
  /**
   * Contains information about the operation used.
   */
  operation: "feature_selection";
  /**
   * Contains information about the specific type of the operation used.
   */
  operationType: "filter_based";
  /**
   * unit input (type to be specified by the child units)
   */
  inputData: {
    /**
     * number of features to select for model training. If equal to 0, will use all available features
     */
    nFeatures: number;
    /**
     * feature selection algorithm following sklearn.feature_selection
     */
    algorithm: "f_regression";
  };
  /**
   * type of the unit
   */
  type: "processing";
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/ml/unit/processing/feature_selection.json */

export type SoftwareDirectoryMlUnitProcessingFeatureSelection = {
  /**
   * Contains information about the operation used.
   */
  operation: "feature_selection";
  /**
   * Contains information about the specific type of the operation used.
   */
  operationType: "filter_based";
  /**
   * unit input (type to be specified by the child units)
   */
  inputData: {
    /**
     * number of features to select for model training. If equal to 0, will use all available features
     */
    nFeatures: number;
    /**
     * feature selection algorithm following sklearn.feature_selection
     */
    algorithm: "f_regression";
  };
  /**
   * type of the unit
   */
  type: "processing";
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
};
 
/** Schema dist/js/schema/software_directory/ml/unit/processing.json */

export type SoftwareDirectoryMlUnitProcessing =
  | {
      /**
       * Contains information about the operation used.
       */
      operation: "data_transformation";
      /**
       * Contains information about the specific type of the operation used.
       */
      operationType: "scale_and_reduce";
      /**
       * unit input (type to be specified by the child units)
       */
      inputData: {
        /**
         * type of scaler to be applied
         */
        scaler: "standard_scaler";
        /**
         * per-feature scaling data
         */
        perFeature?: {
          /**
           * variance in original training data
           */
          variance?: number;
          /**
           * mean value of the original training data
           */
          mean?: number;
          /**
           * scale multiplier for this feature/property
           */
          scale: number;
          /**
           * feature/property name in 'flattened' format
           */
          name: string;
        }[];
      };
      /**
       * type of the unit
       */
      type: "processing";
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * Contains information about the operation used.
       */
      operation: "feature_selection";
      /**
       * Contains information about the specific type of the operation used.
       */
      operationType: "filter_based";
      /**
       * unit input (type to be specified by the child units)
       */
      inputData: {
        /**
         * number of features to select for model training. If equal to 0, will use all available features
         */
        nFeatures: number;
        /**
         * feature selection algorithm following sklearn.feature_selection
         */
        algorithm: "f_regression";
      };
      /**
       * type of the unit
       */
      type: "processing";
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    };
 
/** Schema dist/js/schema/software_directory/modeling/deepmd.json */

export interface DeePMDAppSchema {
  /**
   * entity name
   */
  name?: "deepmd";
  /**
   * Application's short description.
   */
  summary?: "DeePMD is a deep learning package that is based on neural network fitted first-principles data for many-body potential energy representation and molecular dynamics";
  /**
   * Application version. e.g. 5.3.5
   */
  version?: "2.0.2";
  exec?: "dp" | "lmp" | "python";
  /**
   * The short name of the application. e.g. qe
   */
  shortName?: string;
  /**
   * Application build. e.g. VTST
   */
  build?: string;
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * Whether licensing is present
   */
  isLicensed?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/modeling/espresso/arguments.json */

export interface QuantumEspressoArgumentsSchema {
  /**
   * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
   */
  nimage?: number;
  /**
   * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
   */
  npools?: number;
  /**
   * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
   */
  nband?: number;
  /**
   * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
   */
  ntg?: number;
  /**
   * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
   */
  ndiag?: number;
}
 
/** Schema dist/js/schema/software_directory/modeling/espresso.json */

export interface EspressoAppSchema {
  name?: "espresso";
  summary?: "Quantum Espresso";
  version?:
    | "5.2.1"
    | "5.4.0"
    | "6.0.0"
    | "6.3"
    | "6.4.1"
    | "6.5.0"
    | "6.6.0"
    | "6.7.0"
    | "6.8.0"
    | "7.0"
    | "7.2"
    | "7.3";
}
 
/** Schema dist/js/schema/software_directory/modeling/nwchem.json */

export interface NWChem {
  /**
   * entity name
   */
  name?: "NWChem";
  /**
   * Application's short description.
   */
  summary?: "NWChem: a comprehensive and scalable open-source solution for large scale molecular simulations";
  /**
   * Application version. e.g. 5.3.5
   */
  version?: "6.6" | "7.0.2";
  exec?: "nwchem";
  /**
   * The short name of the application. e.g. qe
   */
  shortName?: string;
  /**
   * Application build. e.g. VTST
   */
  build?: string;
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * Whether licensing is present
   */
  isLicensed?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/modeling/unit/execution.json */

export interface ExecutionUnitSchemaForPhysicsBasedSimulationEnginesDefinedUsingEspressoAsExample {
  /**
   * type of the unit
   */
  type: "execution";
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  executable?: {
    /**
     * The name of the executable. e.g. pw.x
     */
    name: string;
    /**
     * _ids of the application this executable belongs to
     */
    applicationId?: string[];
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  flavor?: {
    /**
     * _id of the executable this flavor belongs to
     */
    executableId?: string;
    /**
     * name of the executable this flavor belongs to
     */
    executableName?: string;
    /**
     * name of the application this flavor belongs to
     */
    applicationName?: string;
    input?: {
      templateId?: string;
      templateName?: string;
      /**
       * name of the resulting input file, if different than template name
       */
      name?: string;
    }[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  /**
   * unit input (type to be specified by the application's execution unit)
   */
  input: (
    | {
        /**
         * Input file name. e.g. pw_scf.in
         */
        name: string;
        /**
         * Content of the input file. e.g. &CONTROL    calculation='scf' ...
         */
        content: string;
      }
    | {
        templateId?: string;
        templateName?: string;
        /**
         * name of the resulting input file, if different than template name
         */
        name?: string;
      }
  )[];
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/modeling/vasp.json */

export interface ViennaAbInitoSimulationPackage {
  /**
   * entity name
   */
  name?: "vasp";
  /**
   * Application's short description.
   */
  summary?: "vienna ab-initio simulation package";
  flavor?: "vasp" | "vasp_nscf" | "vasp_bands";
  /**
   * Application version. e.g. 5.3.5
   */
  version?: "5.3.5";
  exec?: "vasp";
  /**
   * The short name of the application. e.g. qe
   */
  shortName?: string;
  /**
   * Application build. e.g. VTST
   */
  build?: string;
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * Whether licensing is present
   */
  isLicensed?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/scripting/jupyter_lab.json */

export interface JupyterLabApplicationSchema {
  /**
   * entity name
   */
  name?: "jupyterLab";
  flavor?: "notebook";
  /**
   * Application's short description.
   */
  summary?: "Jupyter Lab";
  /**
   * Application version. e.g. 5.3.5
   */
  version?: "0.33.12";
  exec?: "jupyter";
  /**
   * The short name of the application. e.g. qe
   */
  shortName?: string;
  /**
   * Application build. e.g. VTST
   */
  build?: string;
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * Whether licensing is present
   */
  isLicensed?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/scripting/python.json */

export interface PythonProgramingLanguageSchema {
  /**
   * entity name
   */
  name?: "python";
  flavor?: "python2" | "python3";
  /**
   * Application's short description.
   */
  summary?: "Python Script";
  /**
   * Application version. e.g. 5.3.5
   */
  version?: "2.7.5" | "3.6.1";
  exec?: "python";
  /**
   * Optional arguments passed to the Python script
   */
  arguments?: string;
  /**
   * Optional environment variables exported before running the Python script
   */
  environment?: {};
  /**
   * Optional Python dependencies, e.g. amqp==1.4.6
   */
  dependencies?: unknown[];
  /**
   * The short name of the application. e.g. qe
   */
  shortName?: string;
  /**
   * Application build. e.g. VTST
   */
  build?: string;
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * Whether licensing is present
   */
  isLicensed?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/scripting/shell.json */

export interface ShellScriptingLanguageSchema {
  /**
   * entity name
   */
  name?: "shell";
  flavor?: "sh" | "bash" | "zsh" | "csh";
  /**
   * Application's short description.
   */
  summary?: "Shell Script";
  /**
   * Application version. e.g. 5.3.5
   */
  version?: "4.2.46";
  exec?: "sh" | "bash" | "zsh" | "csh";
  /**
   * Optional arguments passed to the Shell script
   */
  arguments?: string;
  /**
   * Optional environment variables exported before running the Shell script
   */
  environment?: {};
  /**
   * The short name of the application. e.g. qe
   */
  shortName?: string;
  /**
   * Application build. e.g. VTST
   */
  build?: string;
  /**
   * Whether advanced compute options are present
   */
  hasAdvancedComputeOptions?: boolean;
  /**
   * Whether licensing is present
   */
  isLicensed?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/software_directory/scripting/unit/execution.json */

export interface ExecutionUnitSchemaForScriptingBasedApplications {
  /**
   * type of the unit
   */
  type: "execution";
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  executable?: {
    /**
     * The name of the executable. e.g. pw.x
     */
    name: string;
    /**
     * _ids of the application this executable belongs to
     */
    applicationId?: string[];
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  flavor?: {
    /**
     * _id of the executable this flavor belongs to
     */
    executableId?: string;
    /**
     * name of the executable this flavor belongs to
     */
    executableName?: string;
    /**
     * name of the application this flavor belongs to
     */
    applicationName?: string;
    input?: {
      templateId?: string;
      templateName?: string;
      /**
       * name of the resulting input file, if different than template name
       */
      name?: string;
    }[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  /**
   * unit input (type to be specified by the application's execution unit)
   */
  input: (
    | {
        /**
         * Input file name. e.g. pw_scf.in
         */
        name: string;
        /**
         * Content of the input file. e.g. &CONTROL    calculation='scf' ...
         */
        content: string;
      }
    | {
        templateId?: string;
        templateName?: string;
        /**
         * name of the resulting input file, if different than template name
         */
        name?: string;
      }
  )[];
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/system/_material.json */

export interface MaterialEntityReferenceSchema {
  /**
   * Material class
   */
  cls?: "Material";
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
}
 
/** Schema dist/js/schema/system/_parent_job.json */

export interface ParentJobEntityReferenceSchema {
  /**
   * Job class
   */
  cls?: "Job";
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
}
 
/** Schema dist/js/schema/system/_project.json */

export interface ProjectEntityReferenceSchema {
  /**
   * Project class
   */
  cls?: "Project";
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
}
 
/** Schema dist/js/schema/system/bankable.json */

export interface BankableSchema {
  /**
   * Identity of the corresponding bank entity
   */
  exabyteId?: string;
  /**
   * Hash string which is calculated based on the meaningful fields of the entity. Used to identify equal entities.
   */
  hash?: string;
}
 
/** Schema dist/js/schema/system/consistency_check.json */

/**
 * The output of consistency checks performed on data adhering to JSON schema, but inconsistent with scientific or logical rules, to show problems in UI.
 */
export interface ConsistencyCheck {
  /**
   * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
   */
  key: string;
  /**
   * Name of the consistency check that is performed, which is listed in an enum.
   */
  name: "default" | "atomsTooClose" | "atomsOverlap";
  /**
   * Severity level of the problem, which is used in UI to differentiate.
   */
  severity: "info" | "warning" | "error";
  /**
   * Message generated by the consistency check describing the problem.
   */
  message: string;
}
 
/** Schema dist/js/schema/system/creator.json */

export interface CreatorEntityReferenceSchema {
  /**
   * Creator class
   */
  cls?: "User";
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
}
 
/** Schema dist/js/schema/system/creator_account.json */

export interface CreatorAccountSchema {
  creatorAccount?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
}
 
/** Schema dist/js/schema/system/database_source.json */

/**
 * information about a database source
 */
export interface DatabaseSourceSchema {
  /**
   * ID string for the materials uploaded from a third party source inside the third party source. For materialsproject.org an example ID is mp-32
   */
  id: string | number;
  /**
   * Third party source name, e.g. materials project, 2dmatpedia, ICSD, etc.
   */
  source: string;
  /**
   * Deprecated. To be removed. A flag that is true when material is initially imported from a third party * (as opposed to being independently designed from scratch).
   */
  origin: boolean;
  /**
   * Original response from external source.
   */
  data?: {};
  /**
   * Digital Object Identifier, e.g. 10.1088/0953-8984/25/10/105506
   */
  doi?: string;
  /**
   * The URL of the original record, e.g. https://next-gen.materialsproject.org/materials/mp-48; ToDo: update to use URI type per https://json-schema.org/understanding-json-schema/reference/string#resource-identifiers
   */
  url?: string;
}
 
/** Schema dist/js/schema/system/defaultable.json */

export interface DefaultableEntitySchema {
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
}
 
/** Schema dist/js/schema/system/description.json */

export interface ExtendedBaseEntitySchema {
  /**
   * entity description
   */
  description?: string;
  descriptionObject?: {};
}
 
/** Schema dist/js/schema/system/entity_reference.json */

export interface EntityReferenceSchema {
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity class
   */
  cls?: string;
  /**
   * entity slug
   */
  slug?: string;
}
 
/** Schema dist/js/schema/system/file_source.json */

/**
 * file source with the information inside
 */
export interface FileSourceSchema {
  /**
   * file extension
   */
  extension?: string;
  /**
   * file name without extension
   */
  filename: string;
  /**
   * file content as raw text
   */
  text: string;
  /**
   * MD5 hash based on file content
   */
  hash: string;
}
 
/** Schema dist/js/schema/system/history.json */

export interface HistorySchema {
  history?: {
    id: string;
    revision: number;
  }[];
}
 
/** Schema dist/js/schema/system/iframe_message.json */

/**
 * communication message between iframe and the parent window.
 */
export interface IframeMessageSchema {
  /**
   * The type of the message to distinguish the direction of the message.
   */
  type: "from-iframe-to-host" | "from-host-to-iframe";
  /**
   * The action to be performed upon receiving the message.
   */
  action: "set-data" | "get-data" | "info";
  /**
   * The content of the message with actual data.
   */
  payload: {};
}
 
/** Schema dist/js/schema/system/in_set.json */

export interface SystemInSetSchema {
  inSet?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
    type?: string;
    index?: number;
  }[];
}
 
/** Schema dist/js/schema/system/is_multi_material.json */

export interface IsMultiSchema {
  isMultiMaterial?: boolean;
}
 
/** Schema dist/js/schema/system/is_outdated.json */

export interface IsOutdatedSchema {
  isOutdated?: boolean;
}
 
/** Schema dist/js/schema/system/job_extended.json */

export interface ExtendedJobSchema {
  mode?: string;
  isExternal?: boolean;
  _materials?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  }[];
  _materialsSet?: {
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity class
     */
    cls?: string;
    /**
     * entity slug
     */
    slug?: string;
  };
  purged?: boolean;
  purgedAt?: number;
  dataset?: {};
}
 
/** Schema dist/js/schema/system/message.json */

/**
 * communication message between Rupy and web application.
 */
export interface MessageSchema {
  header: {
    entity: {
      /**
       * job identifier
       */
      _id: string;
      /**
       * entity name.
       */
      name: "job" | "unit";
      /**
       * unit identifier within the workflow
       */
      flowchartId?: string;
      /**
       * source of the message.
       */
      probe?: "monitor" | "postprocessor";
    };
    /**
     * Rupy-Webapp communication schema version.
     */
    version: string;
    /**
     * Timestamp of the message.
     */
    timestamp: number;
  };
  /**
   * Actual payload of the message.
   */
  payload: {};
}
 
/** Schema dist/js/schema/system/metadata.json */

export interface MetadataSchema {
  metadata?: {};
}
 
/** Schema dist/js/schema/system/name.json */

export interface NameEntitySchema {
  /**
   * entity name
   */
  name?: string;
}
 
/** Schema dist/js/schema/system/owner.json */

export interface EntityOwnerReferenceSchema {
  /**
   * Entity owner class
   */
  cls?: "Account";
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
}
 
/** Schema dist/js/schema/system/path.json */

export interface PathSchema {
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
}
 
/** Schema dist/js/schema/system/path_entity.json */

export interface PathEntitySchema {
  /**
   * entity name
   */
  name?: string;
  /**
   * TODO: Use regex once schema draft version has been updated
   */
  path?: string;
}
 
/** Schema dist/js/schema/system/schema_version.json */

export interface SchemaVersion {
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
}
 
/** Schema dist/js/schema/system/scope.json */

export interface ScopeSchema {
  scope?: string;
}
 
/** Schema dist/js/schema/system/set.json */

export interface EntitySetSchema {
  isEntitySet?: boolean;
  entitySetType?: string;
  entityCls?: string;
}
 
/** Schema dist/js/schema/system/sharing.json */

export interface ExtendedSharingSchema {
  sharedCount?: number;
}
 
/** Schema dist/js/schema/system/soft_removable.json */

export interface SoftRemovableEntitySchema {
  /**
   * Timestamp of the moment when entity was removed
   */
  removedAt?: string;
  /**
   * Identifies that entity was removed
   */
  removed?: boolean;
}
 
/** Schema dist/js/schema/system/status.json */

export interface StatusSchema {
  status?: string;
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
}
 
/** Schema dist/js/schema/system/tags.json */

export interface EntityTagsSchema {
  /**
   * entity tags
   */
  tags?: string[];
}
 
/** Schema dist/js/schema/system/timestampable.json */

export interface TimestampableEntitySchema {
  /**
   * entity creation time
   */
  createdAt?: string;
  /**
   * entity last modification time
   */
  updatedAt?: string;
  createdBy?: string;
  updatedBy?: string;
}
 
/** Schema dist/js/schema/system/use_values.json */

export interface UseValuesSchema {
  useValues?: boolean;
}
 
/** Schema dist/js/schema/workflow/base.json */

export interface BaseWorkflowSchema {
  /**
   * Array of characteristic properties calculated by this workflow (TODO: add enums)
   */
  properties?: (string | {})[];
  /**
   * Whether to use the dataset tab in the job designer. Mutually exclusive with using the materials tab.
   */
  isUsingDataset?: boolean;
  /**
   * Array of workflows with the same schema as the current one.
   */
  workflows?: {}[];
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  metadata?: {};
}
 
/** Schema dist/js/schema/workflow/base_flow.json */

export interface BaseFlow {
  /**
   * subworkflow identity
   */
  _id?: string;
  /**
   * Human-readable name of the subworkflow. e.g. Total-energy
   */
  name: string;
  /**
   * Array of characteristic properties calculated by this subworkflow
   */
  properties?: (string | {})[];
  /**
   * compute parameters
   */
  compute?: {
    /**
     * Name of the submission queues: https://docs.mat3ra.com/infrastructure/resource/queues/. Below enums are for Azure, then AWS circa 2022-08, hence the duplication.
     */
    queue:
      | "D"
      | "OR"
      | "OF"
      | "OFplus"
      | "SR"
      | "SF"
      | "SFplus"
      | "GPOF"
      | "GP2OF"
      | "GP4OF"
      | "GPSF"
      | "GP2SF"
      | "GP4SF"
      | "OR4"
      | "OR8"
      | "OR16"
      | "SR4"
      | "SR8"
      | "SR16"
      | "GOF"
      | "G4OF"
      | "G8OF"
      | "GSF"
      | "G4SF"
      | "G8SF";
    /**
     * number of nodes used for the job inside the RMS.
     */
    nodes: number;
    /**
     * number of CPUs used for the job inside the RMS.
     */
    ppn: number;
    /**
     * Wallclock time limit for computing a job. Clock format: 'hh:mm:ss'
     */
    timeLimit: string;
    /**
     * Convention to use when reasoning about time limits
     */
    timeLimitType?: "per single attempt" | "compound";
    /**
     * Job is allowed to restart on termination.
     */
    isRestartable?: boolean;
    /**
     * Email notification for the job: n - never, a - job aborted, b - job begins, e - job ends. Last three could be combined.
     */
    notify?: string;
    /**
     * Email address to notify about job execution.
     */
    email?: string;
    /**
     * Maximum CPU count per node. This parameter is used to let backend job submission infrastructure know that this job is to be charged for the maximum CPU per node instead of the actual ppn. For premium/fast queues where resources are provisioned on-demand and exclusively per user.
     */
    maxCPU?: number;
    /**
     * Optional arguments specific to using application - VASP, Quantum Espresso, etc. Specified elsewhere
     */
    arguments?: {
      /**
       * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
       */
      nimage?: number;
      /**
       * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
       */
      npools?: number;
      /**
       * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
       */
      nband?: number;
      /**
       * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
       */
      ntg?: number;
      /**
       * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
       */
      ndiag?: number;
    };
    /**
     * Cluster where the job is executed. Optional on create. Required on job submission.
     */
    cluster?: {
      /**
       * FQDN of the cluster. e.g. master-1-staging.exabyte.io
       */
      fqdn?: string;
      /**
       * Job's identity in RMS. e.g. 1234.master-1-staging.exabyte.io
       */
      jid?: string;
    };
    /**
     * Computation error. Optional. Appears only if something happens on jobs execution.
     */
    errors?: {
      /**
       * Domain of the error appearance (internal).
       */
      domain?: "rupy" | "alfred" | "celim" | "webapp";
      /**
       * Should be a short, unique, machine-readable error code string. e.g. FileNotFound
       */
      reason?: string;
      /**
       * Human-readable error message. e.g. 'File Not Found: /home/demo/data/project1/job-123/job-config.json'
       */
      message?: string;
      /**
       * Full machine-readable error traceback. e.g. FileNotFound
       */
      traceback?: string;
    }[];
    /**
     * A Python compatible regex to exclude files from upload. e.g. ^.*.txt& excludes all files with .txt suffix
     */
    excludeFilesPattern?: string;
  } | null;
}
 
/** Schema dist/js/schema/workflow/scope.json */

export interface WorkflowScopeSchema {
  global: {
    [k: string]: unknown;
  };
  local: {
    [k: string]: unknown;
  };
}
 
/** Schema dist/js/schema/workflow/subworkflow/unit.json */

export type WorkflowSubworkflowUnitSchema =
  | {
      /**
       * type of the unit
       */
      type: "io";
      subtype: "input" | "output" | "dataFrame";
      source: "api" | "db" | "object_storage";
      input: (
        | {
            /**
             * rest API endpoint
             */
            endpoint: string;
            /**
             * rest API endpoint options
             */
            endpoint_options: {};
            /**
             * the name of the variable in local scope to save the data under
             */
            name?: string;
            [k: string]: unknown;
          }
        | (
            | {
                /**
                 * IDs of item to retrieve from db
                 */
                ids: string[];
                [k: string]: unknown;
              }
            | {
                /**
                 * db collection name
                 */
                collection: string;
                /**
                 * whether the result should be saved as draft
                 */
                draft: boolean;
                [k: string]: unknown;
              }
          )
        | {
            objectData: {
              /**
               * Object storage container for the file
               */
              CONTAINER?: string;
              /**
               * Name of the file inside the object storage bucket
               */
              NAME?: string;
              /**
               * Object storage provider
               */
              PROVIDER?: string;
              /**
               * Region for the object container specified in Container
               */
              REGION?: string;
              /**
               * Size of the file in bytes
               */
              SIZE?: number;
              /**
               * Unix timestamp showing when the file was last modified
               */
              TIMESTAMP?: string;
            };
            /**
             * if a file with the same filename already exists, whether to overwrite the old file
             */
            overwrite?: boolean;
            /**
             * Relative path to the directory that contains the file.
             */
            pathname?: string;
            /**
             * Basename of the file
             */
            basename?: string;
            /**
             * What kind of file this is, e.g. image / text
             */
            filetype?: string;
            [k: string]: unknown;
          }
      )[];
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "reduce";
      /**
       * corresponding map unit flowchart ID
       */
      mapFlowchartId: string;
      /**
       * input information for reduce unit
       */
      input: {
        /**
         * reduce operation, e.g. aggregate
         */
        operation: string;
        /**
         * arguments which are passed to reduce operation function
         */
        arguments: string[];
      }[];
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "condition";
      /**
       * Input information for condition.
       */
      input: {
        /**
         * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
         */
        scope: string;
        /**
         * Name of the input data. e.g. total_energy
         */
        name: string;
      }[];
      /**
       * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
       */
      statement: string;
      /**
       * Flowchart ID reference for `then` part of the condition.
       */
      then: string;
      /**
       * Flowchart ID reference for `else` part of the condition.
       */
      else: string;
      /**
       * Maximum occurrence of the condition, usable for loops.
       */
      maxOccurrences: number;
      /**
       * Throw exception on reaching to maximum occurence.
       */
      throwException?: boolean;
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "assertion";
      /**
       * The statement to be evaluated
       */
      statement: string;
      /**
       * The error message to be displayed if the assertion fails
       */
      errorMessage?: string;
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "execution";
      application: {
        /**
         * The short name of the application. e.g. qe
         */
        shortName?: string;
        /**
         * Application's short description.
         */
        summary?: string;
        /**
         * Application version. e.g. 5.3.5
         */
        version?: string;
        /**
         * Application build. e.g. VTST
         */
        build?: string;
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * Whether licensing is present
         */
        isLicensed?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        [k: string]: unknown;
      };
      executable?: {
        /**
         * The name of the executable. e.g. pw.x
         */
        name: string;
        /**
         * _ids of the application this executable belongs to
         */
        applicationId?: string[];
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      flavor?: {
        /**
         * _id of the executable this flavor belongs to
         */
        executableId?: string;
        /**
         * name of the executable this flavor belongs to
         */
        executableName?: string;
        /**
         * name of the application this flavor belongs to
         */
        applicationName?: string;
        input?: {
          templateId?: string;
          templateName?: string;
          /**
           * name of the resulting input file, if different than template name
           */
          name?: string;
        }[];
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      /**
       * unit input (type to be specified by the application's execution unit)
       */
      input: {
        [k: string]: unknown;
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "assignment";
      /**
       * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
       */
      input?: {
        /**
         * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
         */
        scope: string;
        /**
         * Name of the input data. e.g. total_energy
         */
        name: string;
      }[];
      /**
       * Name of the global variable. e.g. 'x'
       */
      operand: string;
      /**
       * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
       */
      value: string | boolean | number;
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      scope?: string;
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "processing";
      /**
       * Contains information about the operation used.
       */
      operation: string;
      /**
       * Contains information about the specific type of the operation used.
       */
      operationType: string;
      /**
       * unit input (type to be specified by the child units)
       */
      inputData: {
        [k: string]: unknown;
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    };
 
/** Schema dist/js/schema/workflow/subworkflow.json */

export interface Subworkflow {
  /**
   * Contains the Units of the subworkflow
   */
  units: (
    | {
        /**
         * type of the unit
         */
        type: "io";
        subtype: "input" | "output" | "dataFrame";
        source: "api" | "db" | "object_storage";
        input: (
          | {
              /**
               * rest API endpoint
               */
              endpoint: string;
              /**
               * rest API endpoint options
               */
              endpoint_options: {};
              /**
               * the name of the variable in local scope to save the data under
               */
              name?: string;
              [k: string]: unknown;
            }
          | (
              | {
                  /**
                   * IDs of item to retrieve from db
                   */
                  ids: string[];
                  [k: string]: unknown;
                }
              | {
                  /**
                   * db collection name
                   */
                  collection: string;
                  /**
                   * whether the result should be saved as draft
                   */
                  draft: boolean;
                  [k: string]: unknown;
                }
            )
          | {
              objectData: {
                /**
                 * Object storage container for the file
                 */
                CONTAINER?: string;
                /**
                 * Name of the file inside the object storage bucket
                 */
                NAME?: string;
                /**
                 * Object storage provider
                 */
                PROVIDER?: string;
                /**
                 * Region for the object container specified in Container
                 */
                REGION?: string;
                /**
                 * Size of the file in bytes
                 */
                SIZE?: number;
                /**
                 * Unix timestamp showing when the file was last modified
                 */
                TIMESTAMP?: string;
              };
              /**
               * if a file with the same filename already exists, whether to overwrite the old file
               */
              overwrite?: boolean;
              /**
               * Relative path to the directory that contains the file.
               */
              pathname?: string;
              /**
               * Basename of the file
               */
              basename?: string;
              /**
               * What kind of file this is, e.g. image / text
               */
              filetype?: string;
              [k: string]: unknown;
            }
        )[];
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "reduce";
        /**
         * corresponding map unit flowchart ID
         */
        mapFlowchartId: string;
        /**
         * input information for reduce unit
         */
        input: {
          /**
           * reduce operation, e.g. aggregate
           */
          operation: string;
          /**
           * arguments which are passed to reduce operation function
           */
          arguments: string[];
        }[];
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "condition";
        /**
         * Input information for condition.
         */
        input: {
          /**
           * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
           */
          scope: string;
          /**
           * Name of the input data. e.g. total_energy
           */
          name: string;
        }[];
        /**
         * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
         */
        statement: string;
        /**
         * Flowchart ID reference for `then` part of the condition.
         */
        then: string;
        /**
         * Flowchart ID reference for `else` part of the condition.
         */
        else: string;
        /**
         * Maximum occurrence of the condition, usable for loops.
         */
        maxOccurrences: number;
        /**
         * Throw exception on reaching to maximum occurence.
         */
        throwException?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "assertion";
        /**
         * The statement to be evaluated
         */
        statement: string;
        /**
         * The error message to be displayed if the assertion fails
         */
        errorMessage?: string;
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "execution";
        application: {
          /**
           * The short name of the application. e.g. qe
           */
          shortName?: string;
          /**
           * Application's short description.
           */
          summary?: string;
          /**
           * Application version. e.g. 5.3.5
           */
          version?: string;
          /**
           * Application build. e.g. VTST
           */
          build?: string;
          /**
           * Whether advanced compute options are present
           */
          hasAdvancedComputeOptions?: boolean;
          /**
           * Whether licensing is present
           */
          isLicensed?: boolean;
          /**
           * entity identity
           */
          _id?: string;
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * entity name
           */
          name?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          [k: string]: unknown;
        };
        executable?: {
          /**
           * The name of the executable. e.g. pw.x
           */
          name: string;
          /**
           * _ids of the application this executable belongs to
           */
          applicationId?: string[];
          /**
           * Whether advanced compute options are present
           */
          hasAdvancedComputeOptions?: boolean;
          /**
           * entity identity
           */
          _id?: string;
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
        };
        flavor?: {
          /**
           * _id of the executable this flavor belongs to
           */
          executableId?: string;
          /**
           * name of the executable this flavor belongs to
           */
          executableName?: string;
          /**
           * name of the application this flavor belongs to
           */
          applicationName?: string;
          input?: {
            templateId?: string;
            templateName?: string;
            /**
             * name of the resulting input file, if different than template name
             */
            name?: string;
          }[];
          /**
           * entity identity
           */
          _id?: string;
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * entity name
           */
          name?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
        };
        /**
         * unit input (type to be specified by the application's execution unit)
         */
        input: {
          [k: string]: unknown;
        };
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "assignment";
        /**
         * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
         */
        input?: {
          /**
           * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
           */
          scope: string;
          /**
           * Name of the input data. e.g. total_energy
           */
          name: string;
        }[];
        /**
         * Name of the global variable. e.g. 'x'
         */
        operand: string;
        /**
         * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
         */
        value: string | boolean | number;
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        scope?: string;
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "processing";
        /**
         * Contains information about the operation used.
         */
        operation: string;
        /**
         * Contains information about the specific type of the operation used.
         */
        operationType: string;
        /**
         * unit input (type to be specified by the child units)
         */
        inputData: {
          [k: string]: unknown;
        };
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
  )[];
  model: {
    /**
     * general type of the model, eg. `dft`
     */
    type: string;
    /**
     * general subtype of the model, eg. `lda`
     */
    subtype: string;
    method: {
      /**
       * general type of this method, eg. `pseudopotential`
       */
      type: string;
      /**
       * general subtype of this method, eg. `ultra-soft`
       */
      subtype: string;
      /**
       * Object showing the actual possible precision based on theory and implementation
       */
      precision?: {};
      /**
       * additional data specific to method, eg. array of pseudopotentials
       */
      data?: {};
    };
    [k: string]: unknown;
  };
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  /**
   * Defines whether to store the results/properties extracted in this unit to properties collection
   */
  isDraft?: boolean;
  /**
   * subworkflow identity
   */
  _id?: string;
  /**
   * Human-readable name of the subworkflow. e.g. Total-energy
   */
  name: string;
  /**
   * Array of characteristic properties calculated by this subworkflow
   */
  properties?: (string | {})[];
  /**
   * compute parameters
   */
  compute?: {
    /**
     * Name of the submission queues: https://docs.mat3ra.com/infrastructure/resource/queues/. Below enums are for Azure, then AWS circa 2022-08, hence the duplication.
     */
    queue:
      | "D"
      | "OR"
      | "OF"
      | "OFplus"
      | "SR"
      | "SF"
      | "SFplus"
      | "GPOF"
      | "GP2OF"
      | "GP4OF"
      | "GPSF"
      | "GP2SF"
      | "GP4SF"
      | "OR4"
      | "OR8"
      | "OR16"
      | "SR4"
      | "SR8"
      | "SR16"
      | "GOF"
      | "G4OF"
      | "G8OF"
      | "GSF"
      | "G4SF"
      | "G8SF";
    /**
     * number of nodes used for the job inside the RMS.
     */
    nodes: number;
    /**
     * number of CPUs used for the job inside the RMS.
     */
    ppn: number;
    /**
     * Wallclock time limit for computing a job. Clock format: 'hh:mm:ss'
     */
    timeLimit: string;
    /**
     * Convention to use when reasoning about time limits
     */
    timeLimitType?: "per single attempt" | "compound";
    /**
     * Job is allowed to restart on termination.
     */
    isRestartable?: boolean;
    /**
     * Email notification for the job: n - never, a - job aborted, b - job begins, e - job ends. Last three could be combined.
     */
    notify?: string;
    /**
     * Email address to notify about job execution.
     */
    email?: string;
    /**
     * Maximum CPU count per node. This parameter is used to let backend job submission infrastructure know that this job is to be charged for the maximum CPU per node instead of the actual ppn. For premium/fast queues where resources are provisioned on-demand and exclusively per user.
     */
    maxCPU?: number;
    /**
     * Optional arguments specific to using application - VASP, Quantum Espresso, etc. Specified elsewhere
     */
    arguments?: {
      /**
       * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
       */
      nimage?: number;
      /**
       * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
       */
      npools?: number;
      /**
       * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
       */
      nband?: number;
      /**
       * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
       */
      ntg?: number;
      /**
       * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
       */
      ndiag?: number;
    };
    /**
     * Cluster where the job is executed. Optional on create. Required on job submission.
     */
    cluster?: {
      /**
       * FQDN of the cluster. e.g. master-1-staging.exabyte.io
       */
      fqdn?: string;
      /**
       * Job's identity in RMS. e.g. 1234.master-1-staging.exabyte.io
       */
      jid?: string;
    };
    /**
     * Computation error. Optional. Appears only if something happens on jobs execution.
     */
    errors?: {
      /**
       * Domain of the error appearance (internal).
       */
      domain?: "rupy" | "alfred" | "celim" | "webapp";
      /**
       * Should be a short, unique, machine-readable error code string. e.g. FileNotFound
       */
      reason?: string;
      /**
       * Human-readable error message. e.g. 'File Not Found: /home/demo/data/project1/job-123/job-config.json'
       */
      message?: string;
      /**
       * Full machine-readable error traceback. e.g. FileNotFound
       */
      traceback?: string;
    }[];
    /**
     * A Python compatible regex to exclude files from upload. e.g. ^.*.txt& excludes all files with .txt suffix
     */
    excludeFilesPattern?: string;
  } | null;
}
 
/** Schema dist/js/schema/workflow/unit/assertion.json */

export interface AssertionUnitSchema {
  /**
   * type of the unit
   */
  type: "assertion";
  /**
   * The statement to be evaluated
   */
  statement: string;
  /**
   * The error message to be displayed if the assertion fails
   */
  errorMessage?: string;
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/assignment.json */

export interface AssignmentUnitSchema {
  /**
   * type of the unit
   */
  type: "assignment";
  /**
   * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
   */
  input?: {
    /**
     * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
     */
    scope: string;
    /**
     * Name of the input data. e.g. total_energy
     */
    name: string;
  }[];
  /**
   * Name of the global variable. e.g. 'x'
   */
  operand: string;
  /**
   * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
   */
  value: string | boolean | number;
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  scope?: string;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/base.json */

export interface WorkflowBaseUnitSchema {
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * type of the unit
   */
  type: string;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/condition.json */

export interface ConditionUnitSchema {
  /**
   * type of the unit
   */
  type: "condition";
  /**
   * Input information for condition.
   */
  input: {
    /**
     * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
     */
    scope: string;
    /**
     * Name of the input data. e.g. total_energy
     */
    name: string;
  }[];
  /**
   * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
   */
  statement: string;
  /**
   * Flowchart ID reference for `then` part of the condition.
   */
  then: string;
  /**
   * Flowchart ID reference for `else` part of the condition.
   */
  else: string;
  /**
   * Maximum occurrence of the condition, usable for loops.
   */
  maxOccurrences: number;
  /**
   * Throw exception on reaching to maximum occurence.
   */
  throwException?: boolean;
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/execution.json */

export interface ExecutionUnitSchemaBase {
  /**
   * type of the unit
   */
  type: "execution";
  application: {
    /**
     * The short name of the application. e.g. qe
     */
    shortName?: string;
    /**
     * Application's short description.
     */
    summary?: string;
    /**
     * Application version. e.g. 5.3.5
     */
    version?: string;
    /**
     * Application build. e.g. VTST
     */
    build?: string;
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * Whether licensing is present
     */
    isLicensed?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    [k: string]: unknown;
  };
  executable?: {
    /**
     * The name of the executable. e.g. pw.x
     */
    name: string;
    /**
     * _ids of the application this executable belongs to
     */
    applicationId?: string[];
    /**
     * Whether advanced compute options are present
     */
    hasAdvancedComputeOptions?: boolean;
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  flavor?: {
    /**
     * _id of the executable this flavor belongs to
     */
    executableId?: string;
    /**
     * name of the executable this flavor belongs to
     */
    executableName?: string;
    /**
     * name of the application this flavor belongs to
     */
    applicationName?: string;
    input?: {
      templateId?: string;
      templateName?: string;
      /**
       * name of the resulting input file, if different than template name
       */
      name?: string;
    }[];
    /**
     * entity identity
     */
    _id?: string;
    /**
     * entity slug
     */
    slug?: string;
    systemName?: string;
    consistencyChecks?: {
      /**
       * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
       */
      key: string;
      /**
       * Name of the consistency check that is performed, which is listed in an enum.
       */
      name: "default" | "atomsTooClose" | "atomsOverlap";
      /**
       * Severity level of the problem, which is used in UI to differentiate.
       */
      severity: "info" | "warning" | "error";
      /**
       * Message generated by the consistency check describing the problem.
       */
      message: string;
    }[];
    /**
     * entity's schema version. Used to distinct between different schemas.
     */
    schemaVersion?: string;
    /**
     * entity name
     */
    name?: string;
    /**
     * Identifies that entity is defaultable
     */
    isDefault?: boolean;
    /**
     * names of the pre-processors for this calculation
     */
    preProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the post-processors for this calculation
     */
    postProcessors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the monitors for this calculation
     */
    monitors?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
    /**
     * names of the results for this calculation
     */
    results?: (
      | {
          /**
           * The name of this item. e.g. scf_accuracy
           */
          name: string;
        }
      | string
    )[];
  };
  /**
   * unit input (type to be specified by the application's execution unit)
   */
  input: {
    [k: string]: unknown;
  };
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/input/_input.json */

export interface ExecutionUnitInputSchemaForPhysicsBasedSimulationEngines {
  input?: (
    | {
        /**
         * Input file name. e.g. pw_scf.in
         */
        name: string;
        /**
         * Content of the input file. e.g. &CONTROL    calculation='scf' ...
         */
        content: string;
      }
    | {
        templateId?: string;
        templateName?: string;
        /**
         * name of the resulting input file, if different than template name
         */
        name?: string;
      }
  )[];
}
 
/** Schema dist/js/schema/workflow/unit/input/_inputItem.json */

export interface ExecutionUnitInputItemSchemaForPhysicsBasedSimulationEngines {
  /**
   * Input file name. e.g. pw_scf.in
   */
  name: string;
  /**
   * Content of the input file. e.g. &CONTROL    calculation='scf' ...
   */
  content: string;
}
 
/** Schema dist/js/schema/workflow/unit/input/_inputItemId.json */

export interface ExecutionUnitInputIdItemSchemaForPhysicsBasedSimulationEngines {
  templateId?: string;
  templateName?: string;
  /**
   * name of the resulting input file, if different than template name
   */
  name?: string;
}
 
/** Schema dist/js/schema/workflow/unit/input/_inputItemScope.json */

export interface WorkflowUnitInputSchema {
  /**
   * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
   */
  scope: string;
  /**
   * Name of the input data. e.g. total_energy
   */
  name: string;
}
 
/** Schema dist/js/schema/workflow/unit/input/_map_input/values.json */

export interface UnitValuesSchema {
  values?: string;
}
 
/** Schema dist/js/schema/workflow/unit/input/_map_input.json */

export interface UnitMapInputSchema {
  target?: string;
  values?: (number | string | {})[];
  useValues?: boolean;
  scope?: string;
  name?: string;
}
 
/** Schema dist/js/schema/workflow/unit/io/api.json */

export interface DataIORestAPIInputSchema {
  /**
   * rest API endpoint
   */
  endpoint: string;
  /**
   * rest API endpoint options
   */
  endpoint_options: {};
  /**
   * the name of the variable in local scope to save the data under
   */
  name?: string;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/io/db.json */

export type DataIODatabaseInputOutputSchema =
  | {
      /**
       * IDs of item to retrieve from db
       */
      ids: string[];
      [k: string]: unknown;
    }
  | {
      /**
       * db collection name
       */
      collection: string;
      /**
       * whether the result should be saved as draft
       */
      draft: boolean;
      [k: string]: unknown;
    };
 
/** Schema dist/js/schema/workflow/unit/io/object_storage.json */

export interface ObjectStorageIoSchema {
  objectData: {
    /**
     * Object storage container for the file
     */
    CONTAINER?: string;
    /**
     * Name of the file inside the object storage bucket
     */
    NAME?: string;
    /**
     * Object storage provider
     */
    PROVIDER?: string;
    /**
     * Region for the object container specified in Container
     */
    REGION?: string;
    /**
     * Size of the file in bytes
     */
    SIZE?: number;
    /**
     * Unix timestamp showing when the file was last modified
     */
    TIMESTAMP?: string;
  };
  /**
   * if a file with the same filename already exists, whether to overwrite the old file
   */
  overwrite?: boolean;
  /**
   * Relative path to the directory that contains the file.
   */
  pathname?: string;
  /**
   * Basename of the file
   */
  basename?: string;
  /**
   * What kind of file this is, e.g. image / text
   */
  filetype?: string;
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/io.json */

export interface DataIOUnitSchema {
  /**
   * type of the unit
   */
  type: "io";
  subtype: "input" | "output" | "dataFrame";
  source: "api" | "db" | "object_storage";
  input: (
    | {
        /**
         * rest API endpoint
         */
        endpoint: string;
        /**
         * rest API endpoint options
         */
        endpoint_options: {};
        /**
         * the name of the variable in local scope to save the data under
         */
        name?: string;
        [k: string]: unknown;
      }
    | (
        | {
            /**
             * IDs of item to retrieve from db
             */
            ids: string[];
            [k: string]: unknown;
          }
        | {
            /**
             * db collection name
             */
            collection: string;
            /**
             * whether the result should be saved as draft
             */
            draft: boolean;
            [k: string]: unknown;
          }
      )
    | {
        objectData: {
          /**
           * Object storage container for the file
           */
          CONTAINER?: string;
          /**
           * Name of the file inside the object storage bucket
           */
          NAME?: string;
          /**
           * Object storage provider
           */
          PROVIDER?: string;
          /**
           * Region for the object container specified in Container
           */
          REGION?: string;
          /**
           * Size of the file in bytes
           */
          SIZE?: number;
          /**
           * Unix timestamp showing when the file was last modified
           */
          TIMESTAMP?: string;
        };
        /**
         * if a file with the same filename already exists, whether to overwrite the old file
         */
        overwrite?: boolean;
        /**
         * Relative path to the directory that contains the file.
         */
        pathname?: string;
        /**
         * Basename of the file
         */
        basename?: string;
        /**
         * What kind of file this is, e.g. image / text
         */
        filetype?: string;
        [k: string]: unknown;
      }
  )[];
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/map.json */

export interface MapUnitSchema {
  /**
   * type of the unit
   */
  type: "map";
  /**
   * Id of workflow to run inside map
   */
  workflowId: string;
  /**
   * Input information for map.
   */
  input: {
    /**
     * Name of the target variable to substitute using the values below. e.g. K_POINTS
     */
    target: string;
    /**
     * Scope to retrieve `values` from, global or flowchartId. Optional if `values` is given.
     */
    scope?: string;
    /**
     * Name of the variable inside the scope to retrieve `values` from. Optional if `values` is given.
     */
    name?: string;
    /**
     * Sequence of values for the target Jinja variable. Optional if `scope` and `name` are given. This can be used for map-reduce type parallel execution
     */
    values?: (string | number | {})[];
    useValues?: boolean;
  };
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/processing.json */

export interface ProcessingUnitSchema {
  /**
   * type of the unit
   */
  type: "processing";
  /**
   * Contains information about the operation used.
   */
  operation: string;
  /**
   * Contains information about the specific type of the operation used.
   */
  operationType: string;
  /**
   * unit input (type to be specified by the child units)
   */
  inputData: {
    [k: string]: unknown;
  };
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/reduce.json */

export interface ReduceUnitSchema {
  /**
   * type of the unit
   */
  type: "reduce";
  /**
   * corresponding map unit flowchart ID
   */
  mapFlowchartId: string;
  /**
   * input information for reduce unit
   */
  input: {
    /**
     * reduce operation, e.g. aggregate
     */
    operation: string;
    /**
     * arguments which are passed to reduce operation function
     */
    arguments: string[];
  }[];
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit/runtime/_runtime_item_full_object.json */

export interface FullResultSchema {
  /**
   * The name of this item. e.g. 'my_custom_property. <OTHER FIELDS TO BE ADDED>'
   */
  name: string;
}
 
/** Schema dist/js/schema/workflow/unit/runtime/_runtime_item_name_object.json */

export interface NameResultSchema {
  /**
   * The name of this item. e.g. scf_accuracy
   */
  name: string;
}
 
/** Schema dist/js/schema/workflow/unit/runtime/_runtime_item_string.json */

/**
 * name of runtime item in shortened notation
 */
export type RuntimeItemString = string;
 
/** Schema dist/js/schema/workflow/unit/runtime/runtime_item.json */

export type RuntimeItemSchema =
  | {
      /**
       * The name of this item. e.g. scf_accuracy
       */
      name: string;
    }
  | string;
 
/** Schema dist/js/schema/workflow/unit/runtime/runtime_items.json */

export interface RuntimeItemsSchemaPrePostProcessorsMonitorsResults {
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
}
 
/** Schema dist/js/schema/workflow/unit/subworkflow.json */

export interface SubworkflowUnitSchema {
  /**
   * type of the unit
   */
  type: "subworkflow";
  /**
   * entity identity
   */
  _id?: string;
  isDraft?: boolean;
  /**
   * name of the unit. e.g. pw_scf
   */
  name?: string;
  /**
   * Status of the unit.
   */
  status?: "idle" | "active" | "warning" | "error" | "finished";
  /**
   * Whether this unit is the first one to be executed.
   */
  head?: boolean;
  /**
   * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
   */
  flowchartId: string;
  /**
   * Next unit's flowchartId. If empty, the current unit is the last.
   */
  next?: string;
  /**
   * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
   */
  enableRender?: boolean;
  context?: {};
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  /**
   * names of the pre-processors for this calculation
   */
  preProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the post-processors for this calculation
   */
  postProcessors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the monitors for this calculation
   */
  monitors?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * names of the results for this calculation
   */
  results?: (
    | {
        /**
         * The name of this item. e.g. scf_accuracy
         */
        name: string;
      }
    | string
  )[];
  /**
   * entity tags
   */
  tags?: string[];
  statusTrack?: {
    trackedAt: number;
    status: string;
    repetition?: number;
  }[];
  [k: string]: unknown;
}
 
/** Schema dist/js/schema/workflow/unit.json */

export type WorkflowUnitSchema =
  | {
      /**
       * type of the unit
       */
      type: "io";
      subtype: "input" | "output" | "dataFrame";
      source: "api" | "db" | "object_storage";
      input: (
        | {
            /**
             * rest API endpoint
             */
            endpoint: string;
            /**
             * rest API endpoint options
             */
            endpoint_options: {};
            /**
             * the name of the variable in local scope to save the data under
             */
            name?: string;
            [k: string]: unknown;
          }
        | (
            | {
                /**
                 * IDs of item to retrieve from db
                 */
                ids: string[];
                [k: string]: unknown;
              }
            | {
                /**
                 * db collection name
                 */
                collection: string;
                /**
                 * whether the result should be saved as draft
                 */
                draft: boolean;
                [k: string]: unknown;
              }
          )
        | {
            objectData: {
              /**
               * Object storage container for the file
               */
              CONTAINER?: string;
              /**
               * Name of the file inside the object storage bucket
               */
              NAME?: string;
              /**
               * Object storage provider
               */
              PROVIDER?: string;
              /**
               * Region for the object container specified in Container
               */
              REGION?: string;
              /**
               * Size of the file in bytes
               */
              SIZE?: number;
              /**
               * Unix timestamp showing when the file was last modified
               */
              TIMESTAMP?: string;
            };
            /**
             * if a file with the same filename already exists, whether to overwrite the old file
             */
            overwrite?: boolean;
            /**
             * Relative path to the directory that contains the file.
             */
            pathname?: string;
            /**
             * Basename of the file
             */
            basename?: string;
            /**
             * What kind of file this is, e.g. image / text
             */
            filetype?: string;
            [k: string]: unknown;
          }
      )[];
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "reduce";
      /**
       * corresponding map unit flowchart ID
       */
      mapFlowchartId: string;
      /**
       * input information for reduce unit
       */
      input: {
        /**
         * reduce operation, e.g. aggregate
         */
        operation: string;
        /**
         * arguments which are passed to reduce operation function
         */
        arguments: string[];
      }[];
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "condition";
      /**
       * Input information for condition.
       */
      input: {
        /**
         * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
         */
        scope: string;
        /**
         * Name of the input data. e.g. total_energy
         */
        name: string;
      }[];
      /**
       * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
       */
      statement: string;
      /**
       * Flowchart ID reference for `then` part of the condition.
       */
      then: string;
      /**
       * Flowchart ID reference for `else` part of the condition.
       */
      else: string;
      /**
       * Maximum occurrence of the condition, usable for loops.
       */
      maxOccurrences: number;
      /**
       * Throw exception on reaching to maximum occurence.
       */
      throwException?: boolean;
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "assertion";
      /**
       * The statement to be evaluated
       */
      statement: string;
      /**
       * The error message to be displayed if the assertion fails
       */
      errorMessage?: string;
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "execution";
      application: {
        /**
         * The short name of the application. e.g. qe
         */
        shortName?: string;
        /**
         * Application's short description.
         */
        summary?: string;
        /**
         * Application version. e.g. 5.3.5
         */
        version?: string;
        /**
         * Application build. e.g. VTST
         */
        build?: string;
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * Whether licensing is present
         */
        isLicensed?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        [k: string]: unknown;
      };
      executable?: {
        /**
         * The name of the executable. e.g. pw.x
         */
        name: string;
        /**
         * _ids of the application this executable belongs to
         */
        applicationId?: string[];
        /**
         * Whether advanced compute options are present
         */
        hasAdvancedComputeOptions?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      flavor?: {
        /**
         * _id of the executable this flavor belongs to
         */
        executableId?: string;
        /**
         * name of the executable this flavor belongs to
         */
        executableName?: string;
        /**
         * name of the application this flavor belongs to
         */
        applicationName?: string;
        input?: {
          templateId?: string;
          templateName?: string;
          /**
           * name of the resulting input file, if different than template name
           */
          name?: string;
        }[];
        /**
         * entity identity
         */
        _id?: string;
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * entity name
         */
        name?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
      };
      /**
       * unit input (type to be specified by the application's execution unit)
       */
      input: {
        [k: string]: unknown;
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "assignment";
      /**
       * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
       */
      input?: {
        /**
         * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
         */
        scope: string;
        /**
         * Name of the input data. e.g. total_energy
         */
        name: string;
      }[];
      /**
       * Name of the global variable. e.g. 'x'
       */
      operand: string;
      /**
       * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
       */
      value: string | boolean | number;
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      scope?: string;
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "processing";
      /**
       * Contains information about the operation used.
       */
      operation: string;
      /**
       * Contains information about the specific type of the operation used.
       */
      operationType: string;
      /**
       * unit input (type to be specified by the child units)
       */
      inputData: {
        [k: string]: unknown;
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "map";
      /**
       * Id of workflow to run inside map
       */
      workflowId: string;
      /**
       * Input information for map.
       */
      input: {
        /**
         * Name of the target variable to substitute using the values below. e.g. K_POINTS
         */
        target: string;
        /**
         * Scope to retrieve `values` from, global or flowchartId. Optional if `values` is given.
         */
        scope?: string;
        /**
         * Name of the variable inside the scope to retrieve `values` from. Optional if `values` is given.
         */
        name?: string;
        /**
         * Sequence of values for the target Jinja variable. Optional if `scope` and `name` are given. This can be used for map-reduce type parallel execution
         */
        values?: (string | number | {})[];
        useValues?: boolean;
      };
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    }
  | {
      /**
       * type of the unit
       */
      type: "subworkflow";
      /**
       * entity identity
       */
      _id?: string;
      isDraft?: boolean;
      /**
       * name of the unit. e.g. pw_scf
       */
      name?: string;
      /**
       * Status of the unit.
       */
      status?: "idle" | "active" | "warning" | "error" | "finished";
      /**
       * Whether this unit is the first one to be executed.
       */
      head?: boolean;
      /**
       * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
       */
      flowchartId: string;
      /**
       * Next unit's flowchartId. If empty, the current unit is the last.
       */
      next?: string;
      /**
       * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
       */
      enableRender?: boolean;
      context?: {};
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      /**
       * names of the pre-processors for this calculation
       */
      preProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the post-processors for this calculation
       */
      postProcessors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the monitors for this calculation
       */
      monitors?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * names of the results for this calculation
       */
      results?: (
        | {
            /**
             * The name of this item. e.g. scf_accuracy
             */
            name: string;
          }
        | string
      )[];
      /**
       * entity tags
       */
      tags?: string[];
      statusTrack?: {
        trackedAt: number;
        status: string;
        repetition?: number;
      }[];
      [k: string]: unknown;
    };
 
/** Schema dist/js/schema/workflow.json */

export interface WorkflowSchema {
  /**
   * Array of subworkflows. Subworkflow can be an instance of workflow to allow for nesting
   */
  subworkflows: {
    /**
     * Contains the Units of the subworkflow
     */
    units: (
      | {
          /**
           * type of the unit
           */
          type: "io";
          subtype: "input" | "output" | "dataFrame";
          source: "api" | "db" | "object_storage";
          input: (
            | {
                /**
                 * rest API endpoint
                 */
                endpoint: string;
                /**
                 * rest API endpoint options
                 */
                endpoint_options: {};
                /**
                 * the name of the variable in local scope to save the data under
                 */
                name?: string;
                [k: string]: unknown;
              }
            | (
                | {
                    /**
                     * IDs of item to retrieve from db
                     */
                    ids: string[];
                    [k: string]: unknown;
                  }
                | {
                    /**
                     * db collection name
                     */
                    collection: string;
                    /**
                     * whether the result should be saved as draft
                     */
                    draft: boolean;
                    [k: string]: unknown;
                  }
              )
            | {
                objectData: {
                  /**
                   * Object storage container for the file
                   */
                  CONTAINER?: string;
                  /**
                   * Name of the file inside the object storage bucket
                   */
                  NAME?: string;
                  /**
                   * Object storage provider
                   */
                  PROVIDER?: string;
                  /**
                   * Region for the object container specified in Container
                   */
                  REGION?: string;
                  /**
                   * Size of the file in bytes
                   */
                  SIZE?: number;
                  /**
                   * Unix timestamp showing when the file was last modified
                   */
                  TIMESTAMP?: string;
                };
                /**
                 * if a file with the same filename already exists, whether to overwrite the old file
                 */
                overwrite?: boolean;
                /**
                 * Relative path to the directory that contains the file.
                 */
                pathname?: string;
                /**
                 * Basename of the file
                 */
                basename?: string;
                /**
                 * What kind of file this is, e.g. image / text
                 */
                filetype?: string;
                [k: string]: unknown;
              }
          )[];
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "reduce";
          /**
           * corresponding map unit flowchart ID
           */
          mapFlowchartId: string;
          /**
           * input information for reduce unit
           */
          input: {
            /**
             * reduce operation, e.g. aggregate
             */
            operation: string;
            /**
             * arguments which are passed to reduce operation function
             */
            arguments: string[];
          }[];
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "condition";
          /**
           * Input information for condition.
           */
          input: {
            /**
             * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
             */
            scope: string;
            /**
             * Name of the input data. e.g. total_energy
             */
            name: string;
          }[];
          /**
           * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
           */
          statement: string;
          /**
           * Flowchart ID reference for `then` part of the condition.
           */
          then: string;
          /**
           * Flowchart ID reference for `else` part of the condition.
           */
          else: string;
          /**
           * Maximum occurrence of the condition, usable for loops.
           */
          maxOccurrences: number;
          /**
           * Throw exception on reaching to maximum occurence.
           */
          throwException?: boolean;
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "assertion";
          /**
           * The statement to be evaluated
           */
          statement: string;
          /**
           * The error message to be displayed if the assertion fails
           */
          errorMessage?: string;
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "execution";
          application: {
            /**
             * The short name of the application. e.g. qe
             */
            shortName?: string;
            /**
             * Application's short description.
             */
            summary?: string;
            /**
             * Application version. e.g. 5.3.5
             */
            version?: string;
            /**
             * Application build. e.g. VTST
             */
            build?: string;
            /**
             * Whether advanced compute options are present
             */
            hasAdvancedComputeOptions?: boolean;
            /**
             * Whether licensing is present
             */
            isLicensed?: boolean;
            /**
             * entity identity
             */
            _id?: string;
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * entity name
             */
            name?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            [k: string]: unknown;
          };
          executable?: {
            /**
             * The name of the executable. e.g. pw.x
             */
            name: string;
            /**
             * _ids of the application this executable belongs to
             */
            applicationId?: string[];
            /**
             * Whether advanced compute options are present
             */
            hasAdvancedComputeOptions?: boolean;
            /**
             * entity identity
             */
            _id?: string;
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
          };
          flavor?: {
            /**
             * _id of the executable this flavor belongs to
             */
            executableId?: string;
            /**
             * name of the executable this flavor belongs to
             */
            executableName?: string;
            /**
             * name of the application this flavor belongs to
             */
            applicationName?: string;
            input?: {
              templateId?: string;
              templateName?: string;
              /**
               * name of the resulting input file, if different than template name
               */
              name?: string;
            }[];
            /**
             * entity identity
             */
            _id?: string;
            /**
             * entity slug
             */
            slug?: string;
            systemName?: string;
            consistencyChecks?: {
              /**
               * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
               */
              key: string;
              /**
               * Name of the consistency check that is performed, which is listed in an enum.
               */
              name: "default" | "atomsTooClose" | "atomsOverlap";
              /**
               * Severity level of the problem, which is used in UI to differentiate.
               */
              severity: "info" | "warning" | "error";
              /**
               * Message generated by the consistency check describing the problem.
               */
              message: string;
            }[];
            /**
             * entity's schema version. Used to distinct between different schemas.
             */
            schemaVersion?: string;
            /**
             * entity name
             */
            name?: string;
            /**
             * Identifies that entity is defaultable
             */
            isDefault?: boolean;
            /**
             * names of the pre-processors for this calculation
             */
            preProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the post-processors for this calculation
             */
            postProcessors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the monitors for this calculation
             */
            monitors?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
            /**
             * names of the results for this calculation
             */
            results?: (
              | {
                  /**
                   * The name of this item. e.g. scf_accuracy
                   */
                  name: string;
                }
              | string
            )[];
          };
          /**
           * unit input (type to be specified by the application's execution unit)
           */
          input: {
            [k: string]: unknown;
          };
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "assignment";
          /**
           * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
           */
          input?: {
            /**
             * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
             */
            scope: string;
            /**
             * Name of the input data. e.g. total_energy
             */
            name: string;
          }[];
          /**
           * Name of the global variable. e.g. 'x'
           */
          operand: string;
          /**
           * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
           */
          value: string | boolean | number;
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          scope?: string;
          [k: string]: unknown;
        }
      | {
          /**
           * type of the unit
           */
          type: "processing";
          /**
           * Contains information about the operation used.
           */
          operation: string;
          /**
           * Contains information about the specific type of the operation used.
           */
          operationType: string;
          /**
           * unit input (type to be specified by the child units)
           */
          inputData: {
            [k: string]: unknown;
          };
          /**
           * entity identity
           */
          _id?: string;
          isDraft?: boolean;
          /**
           * name of the unit. e.g. pw_scf
           */
          name?: string;
          /**
           * Status of the unit.
           */
          status?: "idle" | "active" | "warning" | "error" | "finished";
          /**
           * Whether this unit is the first one to be executed.
           */
          head?: boolean;
          /**
           * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
           */
          flowchartId: string;
          /**
           * Next unit's flowchartId. If empty, the current unit is the last.
           */
          next?: string;
          /**
           * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
           */
          enableRender?: boolean;
          context?: {};
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * entity tags
           */
          tags?: string[];
          statusTrack?: {
            trackedAt: number;
            status: string;
            repetition?: number;
          }[];
          [k: string]: unknown;
        }
    )[];
    model: {
      /**
       * general type of the model, eg. `dft`
       */
      type: string;
      /**
       * general subtype of the model, eg. `lda`
       */
      subtype: string;
      method: {
        /**
         * general type of this method, eg. `pseudopotential`
         */
        type: string;
        /**
         * general subtype of this method, eg. `ultra-soft`
         */
        subtype: string;
        /**
         * Object showing the actual possible precision based on theory and implementation
         */
        precision?: {};
        /**
         * additional data specific to method, eg. array of pseudopotentials
         */
        data?: {};
      };
      [k: string]: unknown;
    };
    application: {
      /**
       * The short name of the application. e.g. qe
       */
      shortName?: string;
      /**
       * Application's short description.
       */
      summary?: string;
      /**
       * Application version. e.g. 5.3.5
       */
      version?: string;
      /**
       * Application build. e.g. VTST
       */
      build?: string;
      /**
       * Whether advanced compute options are present
       */
      hasAdvancedComputeOptions?: boolean;
      /**
       * Whether licensing is present
       */
      isLicensed?: boolean;
      /**
       * entity identity
       */
      _id?: string;
      /**
       * entity slug
       */
      slug?: string;
      systemName?: string;
      consistencyChecks?: {
        /**
         * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
         */
        key: string;
        /**
         * Name of the consistency check that is performed, which is listed in an enum.
         */
        name: "default" | "atomsTooClose" | "atomsOverlap";
        /**
         * Severity level of the problem, which is used in UI to differentiate.
         */
        severity: "info" | "warning" | "error";
        /**
         * Message generated by the consistency check describing the problem.
         */
        message: string;
      }[];
      /**
       * entity's schema version. Used to distinct between different schemas.
       */
      schemaVersion?: string;
      /**
       * entity name
       */
      name?: string;
      /**
       * Identifies that entity is defaultable
       */
      isDefault?: boolean;
      [k: string]: unknown;
    };
    /**
     * Defines whether to store the results/properties extracted in this unit to properties collection
     */
    isDraft?: boolean;
    /**
     * subworkflow identity
     */
    _id?: string;
    /**
     * Human-readable name of the subworkflow. e.g. Total-energy
     */
    name: string;
    /**
     * Array of characteristic properties calculated by this subworkflow
     */
    properties?: (string | {})[];
    /**
     * compute parameters
     */
    compute?: {
      /**
       * Name of the submission queues: https://docs.mat3ra.com/infrastructure/resource/queues/. Below enums are for Azure, then AWS circa 2022-08, hence the duplication.
       */
      queue:
        | "D"
        | "OR"
        | "OF"
        | "OFplus"
        | "SR"
        | "SF"
        | "SFplus"
        | "GPOF"
        | "GP2OF"
        | "GP4OF"
        | "GPSF"
        | "GP2SF"
        | "GP4SF"
        | "OR4"
        | "OR8"
        | "OR16"
        | "SR4"
        | "SR8"
        | "SR16"
        | "GOF"
        | "G4OF"
        | "G8OF"
        | "GSF"
        | "G4SF"
        | "G8SF";
      /**
       * number of nodes used for the job inside the RMS.
       */
      nodes: number;
      /**
       * number of CPUs used for the job inside the RMS.
       */
      ppn: number;
      /**
       * Wallclock time limit for computing a job. Clock format: 'hh:mm:ss'
       */
      timeLimit: string;
      /**
       * Convention to use when reasoning about time limits
       */
      timeLimitType?: "per single attempt" | "compound";
      /**
       * Job is allowed to restart on termination.
       */
      isRestartable?: boolean;
      /**
       * Email notification for the job: n - never, a - job aborted, b - job begins, e - job ends. Last three could be combined.
       */
      notify?: string;
      /**
       * Email address to notify about job execution.
       */
      email?: string;
      /**
       * Maximum CPU count per node. This parameter is used to let backend job submission infrastructure know that this job is to be charged for the maximum CPU per node instead of the actual ppn. For premium/fast queues where resources are provisioned on-demand and exclusively per user.
       */
      maxCPU?: number;
      /**
       * Optional arguments specific to using application - VASP, Quantum Espresso, etc. Specified elsewhere
       */
      arguments?: {
        /**
         * Processors can be divided into different `images`, each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
         */
        nimage?: number;
        /**
         * Each image can be subpartitioned into `pools`, each taking care of a group of k-points.
         */
        npools?: number;
        /**
         * Each pool is subpartitioned into `band groups`, each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions).
         */
        nband?: number;
        /**
         * In order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to `task` groups so that each group can process several wavefunctions at the same time.
         */
        ntg?: number;
        /**
         * A further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the `linear-algebra group`, a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by n2, where n is an integer; n2 must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations.
         */
        ndiag?: number;
      };
      /**
       * Cluster where the job is executed. Optional on create. Required on job submission.
       */
      cluster?: {
        /**
         * FQDN of the cluster. e.g. master-1-staging.exabyte.io
         */
        fqdn?: string;
        /**
         * Job's identity in RMS. e.g. 1234.master-1-staging.exabyte.io
         */
        jid?: string;
      };
      /**
       * Computation error. Optional. Appears only if something happens on jobs execution.
       */
      errors?: {
        /**
         * Domain of the error appearance (internal).
         */
        domain?: "rupy" | "alfred" | "celim" | "webapp";
        /**
         * Should be a short, unique, machine-readable error code string. e.g. FileNotFound
         */
        reason?: string;
        /**
         * Human-readable error message. e.g. 'File Not Found: /home/demo/data/project1/job-123/job-config.json'
         */
        message?: string;
        /**
         * Full machine-readable error traceback. e.g. FileNotFound
         */
        traceback?: string;
      }[];
      /**
       * A Python compatible regex to exclude files from upload. e.g. ^.*.txt& excludes all files with .txt suffix
       */
      excludeFilesPattern?: string;
    } | null;
  }[];
  /**
   * Contains the Units of the Workflow
   */
  units: (
    | {
        /**
         * type of the unit
         */
        type: "io";
        subtype: "input" | "output" | "dataFrame";
        source: "api" | "db" | "object_storage";
        input: (
          | {
              /**
               * rest API endpoint
               */
              endpoint: string;
              /**
               * rest API endpoint options
               */
              endpoint_options: {};
              /**
               * the name of the variable in local scope to save the data under
               */
              name?: string;
              [k: string]: unknown;
            }
          | (
              | {
                  /**
                   * IDs of item to retrieve from db
                   */
                  ids: string[];
                  [k: string]: unknown;
                }
              | {
                  /**
                   * db collection name
                   */
                  collection: string;
                  /**
                   * whether the result should be saved as draft
                   */
                  draft: boolean;
                  [k: string]: unknown;
                }
            )
          | {
              objectData: {
                /**
                 * Object storage container for the file
                 */
                CONTAINER?: string;
                /**
                 * Name of the file inside the object storage bucket
                 */
                NAME?: string;
                /**
                 * Object storage provider
                 */
                PROVIDER?: string;
                /**
                 * Region for the object container specified in Container
                 */
                REGION?: string;
                /**
                 * Size of the file in bytes
                 */
                SIZE?: number;
                /**
                 * Unix timestamp showing when the file was last modified
                 */
                TIMESTAMP?: string;
              };
              /**
               * if a file with the same filename already exists, whether to overwrite the old file
               */
              overwrite?: boolean;
              /**
               * Relative path to the directory that contains the file.
               */
              pathname?: string;
              /**
               * Basename of the file
               */
              basename?: string;
              /**
               * What kind of file this is, e.g. image / text
               */
              filetype?: string;
              [k: string]: unknown;
            }
        )[];
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "reduce";
        /**
         * corresponding map unit flowchart ID
         */
        mapFlowchartId: string;
        /**
         * input information for reduce unit
         */
        input: {
          /**
           * reduce operation, e.g. aggregate
           */
          operation: string;
          /**
           * arguments which are passed to reduce operation function
           */
          arguments: string[];
        }[];
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "condition";
        /**
         * Input information for condition.
         */
        input: {
          /**
           * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
           */
          scope: string;
          /**
           * Name of the input data. e.g. total_energy
           */
          name: string;
        }[];
        /**
         * Condition statement. e.g. 'abs(x-total_energy) < 1e-5'
         */
        statement: string;
        /**
         * Flowchart ID reference for `then` part of the condition.
         */
        then: string;
        /**
         * Flowchart ID reference for `else` part of the condition.
         */
        else: string;
        /**
         * Maximum occurrence of the condition, usable for loops.
         */
        maxOccurrences: number;
        /**
         * Throw exception on reaching to maximum occurence.
         */
        throwException?: boolean;
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "assertion";
        /**
         * The statement to be evaluated
         */
        statement: string;
        /**
         * The error message to be displayed if the assertion fails
         */
        errorMessage?: string;
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "execution";
        application: {
          /**
           * The short name of the application. e.g. qe
           */
          shortName?: string;
          /**
           * Application's short description.
           */
          summary?: string;
          /**
           * Application version. e.g. 5.3.5
           */
          version?: string;
          /**
           * Application build. e.g. VTST
           */
          build?: string;
          /**
           * Whether advanced compute options are present
           */
          hasAdvancedComputeOptions?: boolean;
          /**
           * Whether licensing is present
           */
          isLicensed?: boolean;
          /**
           * entity identity
           */
          _id?: string;
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * entity name
           */
          name?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          [k: string]: unknown;
        };
        executable?: {
          /**
           * The name of the executable. e.g. pw.x
           */
          name: string;
          /**
           * _ids of the application this executable belongs to
           */
          applicationId?: string[];
          /**
           * Whether advanced compute options are present
           */
          hasAdvancedComputeOptions?: boolean;
          /**
           * entity identity
           */
          _id?: string;
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
        };
        flavor?: {
          /**
           * _id of the executable this flavor belongs to
           */
          executableId?: string;
          /**
           * name of the executable this flavor belongs to
           */
          executableName?: string;
          /**
           * name of the application this flavor belongs to
           */
          applicationName?: string;
          input?: {
            templateId?: string;
            templateName?: string;
            /**
             * name of the resulting input file, if different than template name
             */
            name?: string;
          }[];
          /**
           * entity identity
           */
          _id?: string;
          /**
           * entity slug
           */
          slug?: string;
          systemName?: string;
          consistencyChecks?: {
            /**
             * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
             */
            key: string;
            /**
             * Name of the consistency check that is performed, which is listed in an enum.
             */
            name: "default" | "atomsTooClose" | "atomsOverlap";
            /**
             * Severity level of the problem, which is used in UI to differentiate.
             */
            severity: "info" | "warning" | "error";
            /**
             * Message generated by the consistency check describing the problem.
             */
            message: string;
          }[];
          /**
           * entity's schema version. Used to distinct between different schemas.
           */
          schemaVersion?: string;
          /**
           * entity name
           */
          name?: string;
          /**
           * Identifies that entity is defaultable
           */
          isDefault?: boolean;
          /**
           * names of the pre-processors for this calculation
           */
          preProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the post-processors for this calculation
           */
          postProcessors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the monitors for this calculation
           */
          monitors?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
          /**
           * names of the results for this calculation
           */
          results?: (
            | {
                /**
                 * The name of this item. e.g. scf_accuracy
                 */
                name: string;
              }
            | string
          )[];
        };
        /**
         * unit input (type to be specified by the application's execution unit)
         */
        input: {
          [k: string]: unknown;
        };
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "assignment";
        /**
         * Input information for assignment. if omitted, means that it is an initialization unit, otherwise it is an assignment.
         */
        input?: {
          /**
           * Scope of the variable. e.g. 'global' or 'flowchart_id_2'
           */
          scope: string;
          /**
           * Name of the input data. e.g. total_energy
           */
          name: string;
        }[];
        /**
         * Name of the global variable. e.g. 'x'
         */
        operand: string;
        /**
         * Value of the variable. The value content could be a simple integer, string or a python expression. e.g. '0' (initialization), 'sin(x)+1' (expression)
         */
        value: string | boolean | number;
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        scope?: string;
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "processing";
        /**
         * Contains information about the operation used.
         */
        operation: string;
        /**
         * Contains information about the specific type of the operation used.
         */
        operationType: string;
        /**
         * unit input (type to be specified by the child units)
         */
        inputData: {
          [k: string]: unknown;
        };
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "map";
        /**
         * Id of workflow to run inside map
         */
        workflowId: string;
        /**
         * Input information for map.
         */
        input: {
          /**
           * Name of the target variable to substitute using the values below. e.g. K_POINTS
           */
          target: string;
          /**
           * Scope to retrieve `values` from, global or flowchartId. Optional if `values` is given.
           */
          scope?: string;
          /**
           * Name of the variable inside the scope to retrieve `values` from. Optional if `values` is given.
           */
          name?: string;
          /**
           * Sequence of values for the target Jinja variable. Optional if `scope` and `name` are given. This can be used for map-reduce type parallel execution
           */
          values?: (string | number | {})[];
          useValues?: boolean;
        };
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
    | {
        /**
         * type of the unit
         */
        type: "subworkflow";
        /**
         * entity identity
         */
        _id?: string;
        isDraft?: boolean;
        /**
         * name of the unit. e.g. pw_scf
         */
        name?: string;
        /**
         * Status of the unit.
         */
        status?: "idle" | "active" | "warning" | "error" | "finished";
        /**
         * Whether this unit is the first one to be executed.
         */
        head?: boolean;
        /**
         * Identity of the unit in the workflow. Used to trace the execution flow of the workflow.
         */
        flowchartId: string;
        /**
         * Next unit's flowchartId. If empty, the current unit is the last.
         */
        next?: string;
        /**
         * Whether Rupy should attempt to use Jinja templating to add context variables into the unit
         */
        enableRender?: boolean;
        context?: {};
        /**
         * entity slug
         */
        slug?: string;
        systemName?: string;
        consistencyChecks?: {
          /**
           * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
           */
          key: string;
          /**
           * Name of the consistency check that is performed, which is listed in an enum.
           */
          name: "default" | "atomsTooClose" | "atomsOverlap";
          /**
           * Severity level of the problem, which is used in UI to differentiate.
           */
          severity: "info" | "warning" | "error";
          /**
           * Message generated by the consistency check describing the problem.
           */
          message: string;
        }[];
        /**
         * entity's schema version. Used to distinct between different schemas.
         */
        schemaVersion?: string;
        /**
         * Identifies that entity is defaultable
         */
        isDefault?: boolean;
        /**
         * names of the pre-processors for this calculation
         */
        preProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the post-processors for this calculation
         */
        postProcessors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the monitors for this calculation
         */
        monitors?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * names of the results for this calculation
         */
        results?: (
          | {
              /**
               * The name of this item. e.g. scf_accuracy
               */
              name: string;
            }
          | string
        )[];
        /**
         * entity tags
         */
        tags?: string[];
        statusTrack?: {
          trackedAt: number;
          status: string;
          repetition?: number;
        }[];
        [k: string]: unknown;
      }
  )[];
  /**
   * Array of characteristic properties calculated by this workflow (TODO: add enums)
   */
  properties?: (string | {})[];
  /**
   * Whether to use the dataset tab in the job designer. Mutually exclusive with using the materials tab.
   */
  isUsingDataset?: boolean;
  /**
   * Array of workflows with the same schema as the current one.
   */
  workflows?: {}[];
  /**
   * entity identity
   */
  _id?: string;
  /**
   * entity slug
   */
  slug?: string;
  systemName?: string;
  consistencyChecks?: {
    /**
     * Key of the property of the entity on which the consistency check is performed in Mongo dot notation, e.g. 'basis.coordinates.1'
     */
    key: string;
    /**
     * Name of the consistency check that is performed, which is listed in an enum.
     */
    name: "default" | "atomsTooClose" | "atomsOverlap";
    /**
     * Severity level of the problem, which is used in UI to differentiate.
     */
    severity: "info" | "warning" | "error";
    /**
     * Message generated by the consistency check describing the problem.
     */
    message: string;
  }[];
  /**
   * entity's schema version. Used to distinct between different schemas.
   */
  schemaVersion?: string;
  /**
   * entity name
   */
  name?: string;
  /**
   * Identifies that entity is defaultable
   */
  isDefault?: boolean;
  metadata?: {};
}
 
