/**
 * @file PenMaterials.hh
 * @author: (modified by) Luis Manzanillas
 * @date 2020, Max Planck Institute for Physics
 */


#ifndef PenMaterials_H
#define PenMaterials_H

#include "globals.hh"


class PenMaterials {
public:
  PenMaterials();
  ~PenMaterials();
  
  void Construct();
  
private:
  G4double lightyield;
    
};

#endif
