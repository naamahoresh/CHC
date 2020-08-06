#ifndef MOO_CMA_PARAMS_H
#define MOO_CMA_PARAMS_H

#include <ostream>
#include <FileUtil/Params.h>

//
// My own derived class for managing my configuration files:
//
class MOCMAParams : public Params 
{
public:
  MOCMAParams( int argc, char **argv ) : Params( argc, argv ) {}

  ~MOCMAParams( ) {}

  void readParams() {
    if ( scanFrom(confFile.c_str()) ) 
	std::cout << "Name of the configuration file: " << confFile << std::endl;
    else {
      std::cerr << "No valid configuration file given! (" << confFile  << ")" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void io( std::istream& is, std::ostream& os, FileUtil::iotype type ) {
    FileUtil::io( is, os, "Task"        , Task 	   , 1u		, type);
    FileUtil::io( is, os, "ArchiveSize" , ArchiveSize , 100u	, type);
    FileUtil::io( is, os, "OffSize"     , OffSize     , ArchiveSize, type);
    FileUtil::io( is, os, "Dimension"   , Dimension   , 10u	, type);   
    FileUtil::io( is, os, "Iterations"  , Iterations  , 500u	, type); 
    FileUtil::io( is, os, "Trials"  , Trials      , 10u	, type); 
    FileUtil::io( is, os, "DisplInt"    , DisplInt    , 100u	, type);   
    FileUtil::io( is, os, "ArchiInt"    , ArchiInt    , 250u	, type);   
    FileUtil::io( is, os, "MinInit"     , MinInit     , -1000.	, type);
    FileUtil::io( is, os, "MaxInit"     , MaxInit     , 1000.	, type);
    FileUtil::io( is, os, "MinInit2"    , MinInit2    , 0.	, type);
    FileUtil::io( is, os, "MaxInit2"    , MaxInit2    , 0.	, type);
    FileUtil::io( is, os, "Seed"       , Seed        , 1u	, type);
    FileUtil::io( is, os, "rotBasis"    , rotBasis    , 0u 	, type); 
    FileUtil::io( is, os, "rpcon1"      , rpcon1      , 1.	, type);
    FileUtil::io( is, os, "rpcon2"      , rpcon2      , 2.	, type);
    FileUtil::io( is, os, "resolution"  , resolution  , 100u	, type);
    
    FileUtil::io( is, os, "fileprefix"  , fileprefix  , (std::string) "mocma", type);
    FileUtil::io( is, os, "dir"         , dir         , (std::string) "results", type);
    
    FileUtil::io( is, os, "penfac"      , penfac      , .000001, type);
    FileUtil::io( is, os, "constrHandl" , constrHandl , false	, type);
    FileUtil::io( is, os, "cutOff"      , cutOff      , false	, type);
    FileUtil::io( is, os, "sigma"        , sigma       , (3.*(MaxInit-MinInit)/5.), type);
    FileUtil::io( is, os, "sigmalower"   , sigmalower  , 0.     , type);

    FileUtil::io( is, os, "useSMeasure"  , useSMeasure , true   , type);
  }

  // Method to show the current content of the class variable: 
  void monitor(std::ostream &os = std::cout) const { 
    os << "Task        " << Task 	<< std::endl;
    os << "ArchiveSize " << ArchiveSize	<< std::endl;
    os << "OffSize     " << OffSize	<< std::endl;
    os << "Dimension   " << Dimension	<< std::endl;
    os << "Iterations  " << Iterations	<< std::endl;
    os << "Trials      " << Trials   	<< std::endl;
    os << "DisplInt    " << DisplInt	<< std::endl;		
    os << "ArchiInt    " << ArchiInt	<< std::endl;
    os << "MinInit     " << MinInit	<< std::endl;
    os << "MaxInit     " << MaxInit	<< std::endl;
    os << "MinInit2    " << MinInit2	<< std::endl;
    os << "MaxInit2    " << MaxInit2	<< std::endl;
    os << "Seed        " << Seed	<< std::endl;
    os << "rotBasis    " << rotBasis	<< std::endl;
    os << "rpcon1      " << rpcon1	<< std::endl;
    os << "rpcon2      " << rpcon2	<< std::endl;
    os << "resolution  " << resolution 	<< std::endl;
    
    os << "fileprefix  " << fileprefix	<< std::endl;
    os << "dir         " << dir      	<< std::endl;
    
    os << "penfac      " << penfac	<< std::endl;
    os << "constrHandl " << constrHandl	<< std::endl;
    os << "cutOff      " << cutOff	<< std::endl;
    os << "useSMeasure " << useSMeasure	<< std::endl;
    os << "sigma       " << sigma	<< std::endl;
    os << "sigmalower  " << sigmalower	<< std::endl;
  }
  
  unsigned Seed;
  unsigned ArchiveSize;
  unsigned OffSize;
  unsigned Dimension;
  unsigned Iterations;	
  unsigned Trials;
  unsigned Task;
  unsigned DisplInt;
  unsigned ArchiInt;
  unsigned resolution;
  unsigned rotBasis;
  double   MinInit;
  double   MaxInit;
  double   MinInit2;
  double   MaxInit2;
  double   rpcon1;
  double   rpcon2;
  double   penfac;
  double   sigma;
  double   sigmalower;
  bool	   constrHandl;
  bool	   cutOff;
  bool     useSMeasure;
  std::string   fileprefix;
  std::string   dir;
};


#endif
