# For more information about using CMake with Android Studio, read the
# documentation: https://d.android.com/studio/projects/add-native-code.html

# Sets the minimum version of CMake required to build the native library.

cmake_minimum_required(VERSION 3.4.1)

# Creates and names a library, sets it as either STATIC
# or SHARED, and provides the relative paths to its source code.
# You can define multiple libraries, and CMake builds them for you.
# Gradle automatically packages shared libraries with your APK.

add_library( # Sets the name of the library.
        native-lib

        # Sets the library as a shared library.
        SHARED

        # Provides a relative path to your source file(s).
        native-lib.cpp
        AutoTune.cpp
        IndexBinaryFromFloat.cpp
        IndexIVFPQ.cpp
        MetaIndexes.cpp
        AuxIndexStructures.cpp
        IndexBinaryHNSW.cpp
        IndexIVFSpectralHash.cpp
        native-lib.cpp
        Clustering.cpp
        IndexBinaryIVF.cpp
        IndexLSH.cpp
        OnDiskInvertedLists.cpp
        FaissException.cpp
        Index.cpp
        IndexPQ.cpp
        PolysemousTraining.cpp
        hamming.cpp
        IndexFlat.cpp
        IndexReplicas.cpp
        ProductQuantizer.cpp
        Heap.cpp
        IndexHNSW.cpp
        IndexScalarQuantizer.cpp
        utils.cpp
        HNSW.cpp
        index_io.cpp
        IndexShards.cpp
        utils_simd.cpp
        IndexBinary.cpp
        IndexIVF.cpp
        InvertedLists.cpp
        VectorTransform.cpp
        IndexBinaryFlat.cpp
        IndexIVFFlat.cpp
        IVFlib.cpp
        WorkerThread.cpp
        Extra.cpp
        sgesvd.cpp
        sgeqr2.cpp
        slarft.cpp
        slarfb.cpp
        sorg2r.cpp
        s_cat.c
        f77_aloc.c
        slange.c
        slascl.c
        slaset.c
        slacpy.c
        sgebrd.c
        sorgbr.c
        sbdsqr.c
        sormbr.c
        sorglq.c
        sgelqf.c
        slarfg.c
        slarf.c
        slassq.c
        slartg.c
        slasv2.c
        slas2.c
        slasr.c
        slasq1.c
        sgebd2.c
        sisnan.c
        slabrd.c
        pow_dd.c
        s_copy.c
        scopy.c
        strmm.c
        snrm2.c
        srot.c
        pow_ri.c
        r_sign.c
        strmv.c
        sormlq.c
        sormqr.c
        sgelq2.c
        slapy2.c
        slasq2.c
        slasrt.c
        sorgl2.c
        ilaslr.c
        slaisnan.c
        ilaslc.c
        slaqr3.c
        slarfp.c
        sorm2l.c
        sorm2r.c
        strexc.c
        sormhr.c
        sorml2.c
        slabad.c
        slahqr.c
        slaqr4.c
        slasq3.c
        slasq4.c
        slasq5.c
        slasq6.c
        slanv2.c
        slaqr2.c
        slaqr5.c
        slaexc.c
        sgehrd.c
        r_lg10.c
        slaqr1.c
        slasy2.c
        slarfx.c
        slahr2.c
        sgehd2.c
        saxpy.c
        )

# Searches for a specified prebuilt library and stores the path as a
# variable. Because CMake includes system libraries in the search path by
# default, you only need to specify the name of the public NDK library
# you want to add. CMake verifies that the library exists before
# completing its build.

find_library( # Sets the name of the path variable.
        log-lib

        # Specifies the name of the NDK library that
        # you want CMake to locate.
        log)

add_library( lib_blas SHARED IMPORTED )
set_target_properties(lib_blas PROPERTIES IMPORTED_LOCATION /media/w1/data/Biometric/UareU/Crossmatch/urusdk-android/Android/Samples/Faiss/app/libs/${ANDROID_ABI}/libcoinblas.so)


add_library( lib_lapack SHARED IMPORTED )
set_target_properties(lib_lapack PROPERTIES IMPORTED_LOCATION /media/w1/data/Biometric/UareU/Crossmatch/urusdk-android/Android/Samples/Faiss/app/libs/${ANDROID_ABI}/libcoinlapack.so)


add_library( lib_metis SHARED IMPORTED )
set_target_properties(lib_metis PROPERTIES IMPORTED_LOCATION /media/w1/data/Biometric/UareU/Crossmatch/urusdk-android/Android/Samples/Faiss/app/libs/${ANDROID_ABI}/libcoinmetis.so)

add_library( lib_mumps SHARED IMPORTED )
set_target_properties(lib_mumps PROPERTIES IMPORTED_LOCATION /media/w1/data/Biometric/UareU/Crossmatch/urusdk-android/Android/Samples/Faiss/app/libs/${ANDROID_ABI}/libcoinmumps.so)

add_library( lib_ipopt SHARED IMPORTED )
set_target_properties(lib_ipopt PROPERTIES IMPORTED_LOCATION /media/w1/data/Biometric/UareU/Crossmatch/urusdk-android/Android/Samples/Faiss/app/libs/${ANDROID_ABI}/libipopt.so)

# Specifies libraries CMake should link to your target library. You
# can link multiple libraries, such as libraries you define in this
# build script, prebuilt third-party libraries, or system libraries.

target_link_libraries( # Specifies the target library.
        native-lib

        #external links
        lib_ipopt
        lib_mumps
        lib_metis
        lib_blas
        lib_lapack


        # Links the target library to the log library
        # included in the NDK.
        ${log-lib})