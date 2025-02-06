
if REPEATMASKER_OFFLINE_SETUP:

    localrules: get_repeatmasker_env_name
    rule get_repeatmasker_env_name:
        """This rule will be failing during the first execution
        (when the Conda environment is build for the first time)
        because RepeatMasker tries downloading the default library
        file from an online resource, e.g.,
        https://www.dfam.org/releases/Dfam_3.8/families/Dfam-RepeatMasker.lib.gz
        The wget call potentially has a long timeout - see in bioconda recipe:
        https://github.com/bioconda/bioconda-recipes/blob/master/recipes/repeatmasker/post-link.sh
        """
        output:
            env_name = DIR_PROC.joinpath(
                "05-prepare", "deployment", "repeatmasker_env_name.txt"
            )
        conda:
            DIR_ENVS.joinpath("biotools", "motif_tools.yaml")
        shell:
            "echo ${{CONDA_PREFIX}} > {output}"


    localrules: setup_repeatmasker_default_library
    rule setup_repeatmasker_default_library:
        """
        """
        input:
            conda_env = rules.get_repeatmasker_env_name.output.env_name,
            default_lib = REPEATMASKER_DEFAULT_LIBRARY_FILE
        output:
            rm_default_ok = DIR_PROC.joinpath(
                "05-prepare", "deployment", "repeatmasker_default_lib.ok"
            )
        run:
            import pathlib as pl
            import subprocess as sp

            lib_source = pl.Path(input.default_lib).resolve(strict=True)
            conda_root = open(input.conda_env).read().strip()
            conda_root = pl.Path(conda_root).resolve(strict=True)
            subfolder = "share/RepeatMasker/Libraries"
            target_name = "RepeatMasker.lib"
            lib_target = conda_root.joinpath(subfolder, target_name)

            assert lib_source.suffix == ".gz"

            # NB: use string cmd here b/c of the redirect '>'
            cmd = " ".join(
                ["gzip", "-c", "-d", str(lib_source), ">", str(lib_target)]
            )
            _ = sp.check_call(cmd, shell=True)

            assert lib_target.is_file()

            this_rule = "05-prepare::deployment::repeatmasker::setup_repeatmasker_default_library"
            with open(output.rm_default_ok, "w") as check_file:
                _ = check_file.write(f"{this_rule}\n{get_timestamp()}\n")
        # END OF RUN BLOCK

    REPEATMASKER_SETUP_OK = [
        rules.setup_repeatmasker_default_library.output.rm_default_ok
    ]

else:

    localrules: repeatmasker_online_setup
    rule repeatmasker_online_setup:
        output:
            ok = DIR_PROC.joinpath("05-prepare", "deployment", "repeatmasker_online.ok")
        conda:
            DIR_ENVS.joinpath("biotools", "motif_tools.yaml")
        shell:
            "touch {output.ok}"

    REPEATMASKER_SETUP_OK = [
        rules.repeatmasker_online_setup.output.ok
    ]


if REPEATMASKER_DFAM_ROOT_PART_FILE is not None:

    rule add_dfam_root_partition:
        """The DFAM root partition file is unzipped
        to a size of ~70G (v3.8), so this rule is not local
        and needs time to complete.
        """
        input:
            rm_setup_ok = REPEATMASKER_SETUP_OK,
            conda_env = rules.get_repeatmasker_env_name.output.env_name,
            dfam_root_file = REPEATMASKER_DFAM_ROOT_PART_FILE
        output:
            rm_dfam_root_ok = DIR_PROC.joinpath(
                "05-prepare", "deployment", "repeatmasker_dfam_root.ok"
            )
        resources:
            mem_mb=lambda wildcards, attempt: 2048 * attempt,
            time_hrs=lambda wildcards, attempt: attempt
        run:
            import pathlib as pl
            import subprocess as sp
            import re

            dfam_file_source = pl.Path(input.dfam_root_file).resolve(strict=True)
            conda_root = open(input.conda_env).read().strip()
            conda_root = pl.Path(conda_root).resolve(strict=True)
            subfolder = "share/RepeatMasker/Libraries/famdb"

            extract_dfam_base_name = re.compile(
                "dfam[0-9]+.*full\.0\.h5"
            )
            mobj = extract_dfam_base_name.search(dfam_file_source.name)
            assert mobj is not None, f"Cannot extract base name: {dfam_file_source.name}"
            assert dfam_file_source.suffix == ".gz"

            start, end = mobj.span()
            dfam_file_name = dfam_file_source.name[start:end]

            dfam_file_target = conda_root.joinpath(subfolder, dfam_file_name)

            # NB: use string cmd here b/c of the redirect '>'
            cmd = " ".join(
                ["gzip", "-c", "-d", str(dfam_file_source), ">", str(dfam_file_target)]
            )
            _ = sp.check_call(cmd, shell=True)

            # check if "mini db" file exists and delete
            mini_db = conda_root.joinpath(subfolder, "min_init.0.h5").resolve()
            if mini_db.is_file():
                _ = sp.check_call(["rm", str(mini_db)])

            this_rule = "05-prepare::deployment::repeatmasker::add_dfam_root_partition"
            with open(output.rm_dfam_root_ok, "w") as check_file:
                _ = check_file.write(f"{this_rule}\n{get_timestamp()}\n")
        # END OF RUN BLOCK

    REPEATMASKER_SETUP_OK.append(rules.add_dfam_root_partition.output.rm_dfam_root_ok)


localrules: dump_testseq_file
rule dump_testseq_file:
    output:
        fasta = DIR_LOCAL_REF.joinpath("testseq.fasta")
    run:
        import random as rand
        rand.seed()
        alphabet = list("ACGT")
        start = "".join(rand.choices(alphabet, k=250))
        middle = "AATAA" * 100
        end = "".join(rand.choices(alphabet, k=250))
        sequence = start + middle + end
        with open(output.fasta, "w") as dump:
            _ = dump.write(f">testseq_1k\n{sequence}\n")
    # END OF RUN BLOCK


rule build_repeatmasker_database:
    input:
        fasta = rules.dump_testseq_file.output.fasta,
        rm_setup_ok = list(map(str, REPEATMASKER_SETUP_OK))
    output:
        rm_db_built_ok = DIR_PROC.joinpath(
            "05-prepare", "deployment", "repeatmasker_db_built.ok"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "motif_tools.yaml")
    params:
        species=REPEATMASKER_SPECIES,
        out_dir=lambda wildcards, output: pathlib.Path(
            output.rm_db_built_ok).with_suffix(".wd")
    shell:
        "RepeatMasker -pa 2 -s -dir {params.out_dir} "
        "-species {params.species} {input.fasta} &> {log}"
            " && "
        "touch {output.rm_db_built_ok}"


rule setup_repeatmasker:
    input:
        ok = rules.build_repeatmasker_database.output.rm_db_built_ok
