*** Test Cases ***

Standard MOS reduction
    Setup MOS
    Run MOS reduction   reducedMOS_standard   ACT-CL_J0034.4+0225_P002131N01
    
MOS reduction on selected slits
    Setup MOS
    Run MOS reduction on selected slits   reducedMOS_selected   ACT-CL_J0034.4+0225_P002131N01   SLIT10

MOS reduction using slit file
    Setup MOS
    Make slit file   10   988   1026
    Run MOS reduction using slit file   reducedMOS_slitfile   ACT-CL_J0034.4+0225_P002131N01   manualSlitLoc.txt

*** Keywords ***

Clean up
#     Remove directory        testsCache/reducedMOS_standard              True
#     Remove directory        testsCache/reducedMOS_selected              True
#     Remove directory        testsCache/reducedMOS_slitfile              True
    Remove files            testsCache/*.log
    Remove file             testsCache/manualSlitLoc.txt


*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/RSSMOSPipelineTests.py
Suite Setup                 Setup MOS
Suite Teardown              Clean up


*** Variables ***
