cat directories | parallel --jobs 8 'echo {}; cd {} ; /Users/christopher.illingworth/Documents/GitHub/IVY/run_ivy --input PopSeqs.dat --times Times.in --model I --error 0.207 > Run_ModelIE.out ; cd ../'
