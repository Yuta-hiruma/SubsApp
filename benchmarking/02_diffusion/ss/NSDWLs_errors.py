" This code checks for errors in the input parameters and returns error messages "

def checkerror (Soil_domain,Soil_parameters):
    if Soil_domain.Depth_top<Soil_domain.Depth_bottom : 
         print("~~~~~ ERROR top is smaller than bottom ~~~~~")
    if Soil_domain.Number_o_soils < len(Soil_domain.Depth_o_soils):
         print("~~~~~ ERROR number of soils larger than number of soil depths~~~~~")
    if Soil_domain.Number_o_soils > len(Soil_domain.Depth_o_soils):
         print("~~~~~ ERROR number of soils smaller than number of soil depths~~~~~")   
    if any (i > Soil_domain.Depth_top for i in Soil_domain.Depth_o_soils):
         print("~~~~~ ERROR soil layer above surface ~~~~~")
    if any (i < Soil_domain.Depth_bottom for i in Soil_domain.Depth_o_soils):
         print("~~~~~ ERROR soil layer below surface ~~~~~")
    for  j in Soil_domain.soils:
        if any (i == j for i in Soil_parameters.keys()):
            j
        else:
            print("~~~~~ ERROR unknown soil ~~~~~")
      
    else:
         print("no other boundary input error found")
         
         
" TIC-TOC timing "  
def tic():
    # Homemade version of matlab tic and toc functions
    import time
    global start_tictoc
    startTime_for_tictoc = time.time()


def toc():
    import time
    if 'start_tictoc' in globals():
        print("Elapsed time " + str(time.time() - start_tictoc)
              + " seconds.")
    else:
        print("Toc: start time not set")