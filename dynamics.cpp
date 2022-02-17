#include "mcMd/user/clusterInfo.h"
#include "mcMd/user/mapping.h"

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;
using namespace Util;


int main (){

   std::ifstream commands ("commands");
   /* Format of commands file:
   *  Input Prefix: IPREFIX
   *  Output Prefix: OPREFIX
   *  Start: int (t0)
   *  Increment: int (del_T)
   *  Final: int (tf)
   *  Cutoff unimers: int (cutoff to allocate the cluster to the set of unimers)
   *  Cutoff for preserved micelle: double (cutoff to identify if the micelle 
   *                                        preserved its identity or not)
   *  Cutoff for fission/fusion: double (cutoff to identify which of the micelles 
   *                                     fused/split)                              
   *  Number of histories: int (number of histories to use while analyzing if 
   *                            there are repeated fission or fusion of the 
   *                            same cluster)                                    
   */

   // summary file outputs the all the dynamic processes taking place from the 
   // previous time step -> next time step
 
   std::string l1 ("Input Prefix"); 
   std::string IPre;
   std::string l2 ("Output Prefix");
   std::string OPre;
   std::string l3 ("Start");
   int T0;
   std::string l4 ("Increment");
   int delT;
   std::string l5 ("Final");
   int Tf;
   std::string l6 ("Cutoff unimers");
   int cutoff_U;
   std::string l7 ("Cutoff for preserved micelle");
   double cutoff_P;
   std::string l8 ("Cutoff for fission/fusion");
   double cutoff_F;
   std::string l9 ("Number of histories");
   int Nh;

   // Reading commands file
   std::string line;
   std::string head;
   std::string ParaInput;
   
   while ( std::getline( commands, line ) ) {

      std::stringstream linestream(line);
      std::getline (linestream, head, ':');
      if (head == l1){
         std::getline(linestream, IPre);
         IPre.erase(0,1);
         std::cout<<l1<<": "<<IPre<<std::endl;
      }
      else if (head == l2) {
         std::getline(linestream, OPre);
         OPre.erase(0,1);
         std::cout<<l2<<": "<<OPre<<std::endl;
      } 
      else if (head == l3) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         T0 = stoi (ParaInput);
         std::cout<<l3<<": "<<T0<<std::endl;
      }
      else if (head == l4) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         delT = stoi (ParaInput);
         std::cout<<l4<<": "<<delT<<std::endl;
      }
      else if (head == l5) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         Tf = stoi (ParaInput);
         std::cout<<l5<<": "<<Tf<<std::endl;
      }  
      else if (head == l6) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         cutoff_U = stoi (ParaInput);
         std::cout<<l6<<": "<<cutoff_U<<std::endl;
      }  
      else if (head == l7) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         cutoff_P = stod (ParaInput);
         std::cout<<l7<<": "<<cutoff_P<<std::endl;
      } 
      else if (head == l8) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         cutoff_F = stod (ParaInput);
         std::cout<<l8<<": "<<cutoff_F<<std::endl;
      } 
      else if (head == l9) {
         std::getline(linestream, ParaInput);
         ParaInput.erase(0,1);
         Nh = stoi (ParaInput);
         std::cout<<l9<<": "<<Nh<<std::endl;
      } 
      else {
         std::cout << "Wrong input format"<<std::endl;
      }

   }

   std::cout<<"======================================"<<std::endl;

   /* This vector will keep track of all the dynamic processes
    * taking place in the micelles
    * Index 0: Total chain insertion events
    * Index 1: Total chain expulsion events
    * Index 2: Total stepwise association events
    * Index 3: Total stepwise dissociation events
    * Index 4: Total fission events
    * Index 5: Total fusion events
    */
   std::vector<double > tally;
   // Initializing this vector to 0
   for (int i = 0; i < 7; i++) {
      tally.push_back (0);
   }

   int iSample = 1;

   std::string inClusters (IPre+".clusters");
   std::string inCOMs (IPre+".COMs");
   std::string inMoments (IPre+".momentTensors");
   std::string outClusters (OPre+".clusters");
   std::string outCOMs (OPre+".COMs");
   std::string outMoments (OPre+".momentTensors");

   // Setting the file streams and cluster variables
   std::ifstream inFileClusters0 (inClusters+to_string(T0));
   std::ifstream inFileCOMs0 (inCOMs+to_string(T0));
   std::ifstream inFileMoments0 (inMoments+to_string(T0));
   std::ifstream inFileClusters1;
   std::ifstream inFileCOMs1;
   std::ifstream inFileMoments1;


   std::ofstream outFileClusters0 (outClusters+to_string(T0));
   std::ofstream outFileCOMs0 (outCOMs+to_string(T0));
   std::ofstream outFileMoments0 (outMoments+to_string(T0));
   std::ofstream outFileClusters1;
   std::ofstream outFileCOMs1;
   std::ofstream outFileMoments1;
   //std::ofstream outFileName0 (OPre+to_string(T0));
   //std::ofstream outFileName1 (OPre+to_string(T0+delT));

   std::ofstream summary ("summary");
   clusterInfo step0;
   clusterInfo step1;
   clusterInfo* add0;
   clusterInfo* add1;

   // Associating pointers with the objects
   //add0 = &step0;
   //add1 = &step1; 

   /* Read the first two time steps using readStep0 so as to allocate 
   * the DArrays
   */
   //step0.readStep0(inFileName0, cutoff_U);
   //step1.readStep0(inFileName1, cutoff_U);

   step0.readStep(inFileClusters0, inFileCOMs0, inFileMoments0, cutoff_U);
   //step1.readStep(inFileClusters1, inFileCOMs1, inFileMoments1, cutoff_U);

   // Update maxClusterId. This is a static variable.
   // After this it needs to be updated in the mapping function
   step0.maxClusterId = step0.nClusters;

   std::cout<<"Outside loop ----- Max Cluster ID :"<<step0.maxClusterId<<std::endl;

   // Checking if nClusters at step1 is greater than maxId set till now
   //maxId (& step0, & step1);

   //summary<<"SAMPLE : "<<iSample<<std::endl<<std::endl;

   // Add a MAPPING command over here
   //mapping (& step0, & step1, cutoff_P, cutoff_F, & tally, summary);   
 
   // Writing these timesteps to new output files
   step0.writeStep(outFileClusters0, outFileCOMs0, outFileMoments0);
   //step1.writeStep(outFileName1);
   
   //summary<<std::endl<<std::endl; 

   // Clearing the input and output streams for step0. Reinitialize 
   // using the open function of these streams.
   //inFileClusters0.close(); 
   //inFileCOMs0.close();
   //inFileMoments0.close();

   //outFileClusters0.close();
   //outFileCOMs0.close();
   //outFileMoments0.close();

   // iFile is the string that needs to be appended to the filename
   // iSample is the number of samples that have been analyzed
   // 2 files have already been mapped using readStep0
   //iSample++;

   for (int iFile = (T0+delT); iFile <= Tf; iFile+=delT) {

      // Reinitializing the IO variable

      inFileClusters1.open (inClusters+to_string(iFile));
      inFileCOMs1.open (inCOMs+to_string(iFile));
      inFileMoments1.open (inMoments+to_string(iFile));

      outFileClusters1.open (outClusters+to_string(iFile));
      outFileCOMs1.open (outCOMs+to_string(iFile));
      outFileMoments1.open (outMoments+to_string(iFile));

      //inFileName1.open(IPre+to_string(iFile));
      //outFileName1.open(OPre+to_string(iFile));
      summary<<"SAMPLE : "<<iSample<<std::endl;
      summary<<std::endl;

      // Reading the next cluster file
      step1.readStep(inFileClusters1, inFileCOMs1, inFileMoments1, cutoff_U);

      std::cout<<"************iFile : "<<iFile<<"*************"<<std::endl; 

      // Checking if nClusters at step1 is greater than maxId set till now
      // Why is it required ?
      // => In the function mapping, whenever there is a birth/death of a new cluster
      // maxClusterId is incremented and the corresponding cluster at step1 is renamed.
      // Now assume that maxID was not called and nClusters at step0 is 4 and nClusters
      // at step 1 is 7. So in this case a new cluster at step1 might be given the name
      // "5" which would clash with the cluster name "5" which would already exist at 
      // step1. Therefore, the function clusterIndex in clusterInfo.cpp would give the 
      // wrong index if you call clusterIndex (5).
      maxId (& step0, & step1);

      // MAPPING from step0 to step1
      mapping (& step0, & step1, cutoff_P, cutoff_F, & tally, summary);

      // Writing the step1 file 
      step1.writeStep(outFileClusters1, outFileCOMs1, outFileMoments1);
      summary<<std::endl<<std::endl;

      add1 = &step1;
      add0 = &step0;

      iSample++;

      // Clearing the IO variables for step 1
      inFileClusters1.close();
      inFileCOMs1.close();
      inFileMoments1.close();

      outFileClusters1.close();
      outFileCOMs1.close();
      outFileMoments1.close();

      // Copying step1 to step0 after clearing step0.
      step0 = *add1;
      step1 = *add0;
      step1.clear();
   }

   std::cout<<"======================================"<<std::endl;
   std::cout<<"Total samples analyzed: " << iSample << std::endl;
 
   std::cout<<std::endl<<"Thank you.....Have a nice day....:D"<<std::endl;
 
   summary<<"============================================="<<std::endl<<std::endl;
   summary<<"Total samples analyzed: " << iSample << std::endl;

   for (int i = 0; i < tally.size()-1; i++) {
      tally [i] = tally[i]/((double)iSample-1.0);
   } 
   summary<<"Chain insertion events per step : "<< tally [0] << setprecision(8) << std::endl;
   summary<<"Chain expulsion events per step : "<< tally [1] << setprecision(8) << std::endl;
   summary<<"Stepwise association events per step : "<< tally [2] << setprecision(8) << std::endl;
   summary<<"Stepwise dissociation events per step : "<< tally [3] << setprecision(8) << std::endl;
   summary<<"Fission events per step : "<< tally [4] << setprecision(8) << std::endl;
   summary<<"Fusion events per step : "<< tally [5] << setprecision(8) << std::endl;
   summary<<"Total algorithmic errors : "<< tally [6] << setprecision(8) << std::endl;
 
   //std::cout<<"Is this it....?"<< std::endl;


   // Do the initial allocation of maxClusterId using the nClusters after
   // calling readStep0 for step0
   // Start incrementing it in mapping thereafter

   // Take the initial timestep, final timestep and the increment as an input

   // Then have a for loop which reads all the files one by one, maps the 
   // clusters and then writes it to a seperate folder. You might also need 
   // to have the path of this folder as an input. 

   // need to have a function in this file to map the clusters
 
   // might also need to write a function to switch step0 and step1 as 
   // when you simply copy, the inequality between the sizes of the two 
   // GArrays/ DArrays might turn up as an error.

   return 0;
}
