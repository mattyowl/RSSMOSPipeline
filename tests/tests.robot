*** Test Cases ***

Standard longslit reduction
    Setup longslit
    Run reduction on selected slits      reduced_long_slit   all    SLIT7

No flats longslit reduction
    Setup no flats longslit
    Run no flats reduction on selected slits      reduced_long_slit_no_flats   all    SLIT7

Standard MOS reduction
    Setup MOS
    Run reduction   reduced_MOS_standard   ACT-CL_J0034.4+0225_P002131N01
    
MOS reduction on selected slits
    Setup MOS
    Run reduction on selected slits   reduced_MOS_selected   ACT-CL_J0034.4+0225_P002131N01   SLIT10

MOS reduction using slit file
    Setup MOS
    Make slit file   10   988   1026
    Run reduction using slit file   reduced_MOS_slitfile   ACT-CL_J0034.4+0225_P002131N01   manualSlitLoc.txt

*** Keywords ***

Clean up
#     Remove directory        testsCache/reduced_MOS_standard              True
#     Remove directory        testsCache/reduced_MOS_selected              True
#     Remove directory        testsCache/reduced_MOS_slitfile              True
    Remove files            testsCache/*.log
    Remove file             testsCache/manualSlitLoc.txt


*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/RSSMOSPipelineTests.py
Suite Setup                 Setup MOS
Suite Teardown              Clean up


*** Variables ***
