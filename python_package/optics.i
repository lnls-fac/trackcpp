class Twiss {

public:
  Pos<double> co;
  Vector etax;
  Vector etay;
  double mux, betax, alphax;
  double muy, betay, alphay;
  Twiss() : co(std::nan("")),
            etax(Vector({std::nan(""),std::nan("")})),
            etay(Vector({std::nan(""),std::nan("")})),
            mux(0), betax(std::nan("")), alphax(std::nan("")),
            muy(0), betay(std::nan("")), alphay(std::nan("")) {}

  bool isundef() const { return std::isnan(this->betax); }

};

Status::type calc_twiss(const Accelerator& accelerator,
                        const Pos<double>& fixed_point,
                        Matrix& m66,
                        std::vector<Twiss>& twiss,
                        Twiss twiss0 = Twiss(),
                        bool closed_flag = false);
