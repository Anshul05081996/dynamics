#ifndef READ_SAMPLE_H
#define READ_SAMPLE_H

#include <vector>

class readSample {

   public:

      // Will have to clear these vectors as well when 
      // appending the vectors 
      void clear();
 
      // The size of this vector will be nFis by 4
      // where nFis will be the number of fission 
      // events for this sample.
      // The first column will have the cluster Id
      // undergoing fision. The next two columns will
      // have the clusterIds being formed by the 
      // cluster undergoing fission. 
      // The third column will be used as a boolean coloumn. 
      // If it has 0, that means the corresponding fission 
      // process has not been reverted.
      std::vector< std::vector<int > > fission;

      // The size of this vector will be nFus by 4
      // where nFus will be the number of fusion 
      // events for this sample.
      // The first two columns will have the cluster Id
      // undergoing fusion. The next column will have
      // the clusterId being formed by the fusion 
      // process. 
      // The third column will be used as a boolean coloumn. 
      // If it has 0, that means the corresponding fission 
      // process has not been reverted.
      // The fourth column will be used as a boolean coloumn.
      // It will be used to check if the chains mixed after
      // a fusion and a corresponding fission. This will be 
      // done by comparing the molecule ids. If there was no
      // mixing the files will be rewritten by bringing back
      // the old cluster ids. If there was mixing, This
      // process will be reported as mixing fusion and the 
      // corresponding fission as mixing fission
      std::vector< std::vector<int > > fusion;   

      // This is the sample number being processed
      int iSample;

};

inline 
void readSample::clear()
{   
   for (int i = 0; i < fission.size(); i++) {
      fission [i].clear();
   }  
   fission.clear();

   for (int i = 0; i < fusion.size(); i++) {
      fusion [i].clear();
   }   
   fusion.clear();

   iSample = -1;
}  

#endif
