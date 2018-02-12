#include <iostream>
#include <string>

// llcv
#include "LLCVBase/Processor.h"
#include "InterTool_App/InterModule.h"

int main( int nargs, char** argv ) {

  std::cout << "[ Select Tune Sample ]" << std::endl;

  std::string cfg        = argv[1];
  std::string INTER_FILE = argv[2];
  std::string SSNET_FILE = argv[3];
  std::string VTX_FILE   = argv[4];
  std::string SHR_FILE   = argv[5];
  std::string TRK_FILE   = argv[6];
  std::string FLSH_FILE  = argv[7];
  
  llcv::Processor llcv_proc;
  llcv::InterModule llcvmod;

  llcv::InterDriver& driver = llcvmod.Driver();
  driver.AttachInterFile(INTER_FILE,"vertex_tree");
  driver.SetOutputFilename("aho_interfile.root");
  //driver.AddSelection(llcv.InterSelToy());

  llcv_proc.add_llcv_ana(&llcvmod);

  llcv_proc.configure( cfg );

  llcv_proc.set_output_lcv_name( "aho_larcv.root" );
  llcv_proc.set_output_ll_name(  "aho_larlite.root" );

  llcv_proc.add_lcv_input_file(SSNET_FILE);
  llcv_proc.add_lcv_input_file(VTX_FILE);
  llcv_proc.add_ll_input_file(SHR_FILE);
  llcv_proc.add_ll_input_file(TRK_FILE);
  llcv_proc.add_ll_input_file(FLSH_FILE);

  llcv_proc.initialize();

  //proc.batch_process_lcv_reverse()

  llcv_proc.finalize();

  return 0;
}
