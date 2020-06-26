workflow test_location {
    call find_tools
}

task find_tools {
    command {
        ls $STAR_ROOT
        echo "@@@@@@@@@@@@@@@@"
        ls $PICARD_ROOT
        echo "@@@@@@@@@@@@@@@@"
        ls $JAVA_ROOT
        echo "@@@@@@@@@@@@@@@@"

        echo $PATH
        echo "################"
        echo $MANPATH
        echo "################"
        echo $LD_LIBRARY_PATH
        echo "################"
    }
    output{
        String message = read_string(stdout())
    }
    runtime {
        docker: "g3chen/starmulti@sha256:0c35fb8426d1045d55544aec97305ec584ec3dfb80d41ccf1adba39a773bec04"
        modules: "star/2.7.3a picard/2.19.2"
    }
}
