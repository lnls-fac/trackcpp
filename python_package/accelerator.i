class Accelerator {
public:
  Accelerator(const double& energy=-1);        // energy < electron_rest_energy -> energy = electron_rest_energy
  double                  energy;              // [eV]
  bool                    cavity_on;
  bool                    radiation_on;
  bool                    vchamber_on;
  int                     harmonic_number;
  std::vector<Element>    lattice;
  std::vector<Kicktable*> kicktables;

  bool operator==(const Accelerator& o) const;
  bool operator!=(const Accelerator& o) const { return !(*this == o); };
  bool isequal(const Accelerator& a) const { return *this == a; }

};
