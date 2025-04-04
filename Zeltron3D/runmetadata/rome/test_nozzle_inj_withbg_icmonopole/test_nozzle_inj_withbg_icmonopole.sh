SIM_ROOT=${CCCSCRATCHDIR}/jetsims_3d
SIM=test_nozzle_inj_withbg_icmonopole
BASE_DIR=$(dirname ${BASH_SOURCE[0]})
RUN_CMD="ccc_msub ./submit_rome.sh"

# Define any extra files that need to be copied to the test simulation
# directory (EXTRA_MOVES), as well as what they should be named once they get
# there (EXTRA_MOVES_TO).
MOVES=( ${SIM}_input.f90 ${SIM}_Makefile ${SIM}_run.sh )
MOVES_TO=( mod_input.f90 Makefile submit_rome.sh )

if [ "$#" -ne 2 ]; then
    echo "Script requires 2 arguments."
    echo "usage: ${0} <zeltron_source_dir> [ <scripmode (one of {'build', 'run'})> ]"
    exit 0
fi

SOURCE_DIR=$1

if [ "$2" == "" ] || [ "$2" == "build" ]; then
    if [ ${#EXTRA_MOVES[@]} != ${#EXTRA_MOVES_TO[@]} ]; then
	echo "${0}: Check my code to make sure EXTRA_MOVES and EXTRA_MOVES_TO"
	echo "have the same length."
	exit
    fi

    if [ ! -d ${SIM_ROOT} ]; then
        echo "Directory ${SIM_ROOT} not found. Please create it first."
        exit
    fi

    # In build mode, we create and populate the base-level simulation directory
    # where the run will take place. This population involves two steps:
    # 1) We copy ALL the Zeltron source-code files from the repository to the
    #    directory where the run will occur.
    # 2) We overwite the few source files whose contents need to be modified
    #    from those of the main/central repository in order to achieve the
    #    desired run. These are the files specified by the MOVES and MOVES_TO
    #    variables.
	RUN_DIR=${SIM_ROOT}/${SIM}
    echo ${RUN_DIR}

	echo "Creating ${RUN_DIR}"
	mkdir ${RUN_DIR}

    # Do step 1) described above.
	echo "Populating ${RUN_DIR}"
	cp ${SOURCE_DIR}/* ${RUN_DIR}
    # mkdir ${RUN_DIR}/2d
    # cp ${SOURCE_DIR}/2d/* ${RUN_DIR}/2d/

    cd ${RUN_DIR}
    bash ./createDataDirs
    bash ./clean
    cd -

    # Now do step 2) described above.
	for i in "${!MOVES[@]}"; do
	    cp ${BASE_DIR}/${MOVES[$i]} ${RUN_DIR}/${MOVES_TO[$i]}
	done

	echo "Building ${RUN_DIR}"
	cd ${RUN_DIR}
    # Compile the code
	# bash ./script_compile
    make
    cd -
elif [ "$2" == "run" ]; then
	RUN_DIR=${SIM_ROOT}/${SIM}

    # Invoking this script in 'run' mode assumes it has already been run in
    # 'build' mode. This mode simply starts the simulation in the appropriate
    # run directory.
	echo "Running zeltron.exe in ${RUN_DIR}"
	cd ${RUN_DIR}
    echo "Current working directory is now:"
    pwd
    echo "Executing command: ${RUN_CMD}"
    # By default, output of the simulation is redirected to zeltron.out and
    # zeltron.err. If you prefer for Zeltron output to be sent to the terminal,
    # you can comment out the following line.
	exec > zeltron.out 2> zeltron.err
    $RUN_CMD

	cd ${BASE_DIR}
else
    echo "Argument '${2}' not recognized"
    echo "usage: ${0} <zeltron_source_dir> [ <scripmode (one of {'build', 'run'})> ]"
fi
