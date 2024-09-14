TEMPLATE = app
CONFIG += console c++20
CONFIG -= app_bundle
CONFIG -= qt

CONFIG(release, debug|release) {
    DEFINES += NDEBUG # makes assert(...) do nothing
}

# DEFINES += DEBUG_VERBOSE
# DEFINES += DEBUG_VERBOSE_ENC_DEC

# choose one from the next two for systematic codewords

DEFINES += PACK_UNSCRAMBLED_PARITY_FIRST
# DEFINES += PACK_UNSCRAMBLED_DATA_FIRST
# DEFINES += PACK_SCRAMBLED // TODO: non-systematic codewords

# TODO: move packing strategy from compiler options to the user code

# DEFINES += USE_BRUTE_FORCE_ELP_ROOTS_SEARCH # default is Chien's Search

SOURCES += \
        main.cpp \
        tests.cpp

HEADERS += \
    bch.h \
    galois.h \
    gf2_polynomial.h \
    measure_time.h \
    polynomial.h \
    tests.h

DISTFILES += \
    LICENSE \
    README.md
