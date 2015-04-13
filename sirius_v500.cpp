#include "lattice.h"
#include "elements.h"
#include <algorithm>
#include <cmath>

void sirius_v500(std::vector<Element>& the_ring) {

	int  harmonic_number = 864;

	//double energy = 3e9;

	// AC10_5
	double qaf_strength       =  2.536876;
	double qad_strength       = -2.730416;
	double qbd2_strength      = -3.961194;
	double qbf_strength       =  3.902838;
	double qbd1_strength      = -2.966239;
	double qf1_strength       =  2.367821;
	double qf2_strength       =  3.354286;
	double qf3_strength       =  3.080632;
	double qf4_strength       =  2.707639;
	double sa1_strength       = -115.7829759411277/2;
	double sa2_strength       =   49.50386128829739/2;
	double sb1_strength       = -214.5386552515188/2;
	double sb2_strength       =  133.1252391065637/2;
	double sd1_strength       = -302.6188062085843/2;
	double sf1_strength       =  369.5045185071228/2;
	double sd2_strength       = -164.3042864671946/2;
	double sd3_strength       = -289.9270429064217/2;
	double sf2_strength       =  333.7039740852999/2;


    //""" --- drift spaces --- """

    double id_length = 2.0; // [m]

    Element dia1     = Element::drift("dia", id_length/2);
    Element dia2     = Element::drift("dia", 3.26920 + 3.65e-3 - id_length/2);
    Element dib1     = Element::drift("dib", id_length/2);
    Element dib2     = Element::drift("dib", 2.909200 + 3.65e-3 - id_length/2);
    Element d10      = Element::drift("d10", 0.100000);
    Element d11      = Element::drift("d11", 0.110000);
    Element d12      = Element::drift("d12", 0.120000);
    Element d13      = Element::drift("d13", 0.130000);
    Element d15      = Element::drift("d15", 0.150000);
    Element d17      = Element::drift("d17", 0.170000);
    Element d18      = Element::drift("d18", 0.180000);
    Element d20      = Element::drift("d20", 0.200000);
    Element d22      = Element::drift("d22", 0.220000);
    Element d23      = Element::drift("d23", 0.230000);
    Element d26      = Element::drift("d26", 0.260000);
    Element d32      = Element::drift("d32", 0.320000);
    Element d44      = Element::drift("d44", 0.440000);

    //""" --- markers --- """
    Element mc     = Element::marker("mc");
	Element mia    = Element::marker("mia");
	Element mib    = Element::marker("mib");
	Element mb1    = Element::marker("mb1");
	Element mb2    = Element::marker("mb2");
	Element mb3    = Element::marker("mb3");
	Element inicio = Element::marker("inicio");
	Element fim    = Element::marker("fim");
	Element mida   = Element::marker("id_enda");
	Element midb   = Element::marker("id_endb");

    //""" --- beam position monitors --- """
	Element mon    = Element::marker("BPM");

    //""" --- quadrupoles --- """
    Element qaf      = Element::quadrupole("qaf",  0.340000, qaf_strength);
	Element qad      = Element::quadrupole("qad",  0.140000, qad_strength);
	Element qbd2     = Element::quadrupole("qbd2", 0.140000, qbd2_strength);
	Element qbf      = Element::quadrupole("qbf",  0.340000, qbf_strength);
	Element qbd1     = Element::quadrupole("qbd1", 0.140000, qbd1_strength);
	Element qf1      = Element::quadrupole("qf1",  0.250000, qf1_strength);
	Element qf2      = Element::quadrupole("qf2",  0.250000, qf2_strength);
	Element qf3      = Element::quadrupole("qf3",  0.250000, qf3_strength);
	Element qf4      = Element::quadrupole("qf4",  0.250000, qf4_strength);

    //""" --- bending magnets --- """

    double deg_2_rad = (M_PI/180.0);

    std::string dip_nam;
	double      dip_len, dip_ang, dip_K, dip_S;
    Element     h1, h2;

    //""" -- b1 -- """
    dip_nam =  "b1";
    dip_len =  0.828080;
    dip_ang =  2.766540 * deg_2_rad;
    dip_K   = -0.78;
    dip_S   =  0;
    h1      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 1*dip_ang/2, 0*dip_ang/2, dip_K, dip_S);
    h2      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 0*dip_ang/2, 1*dip_ang/2, dip_K, dip_S);
    std::vector<Element> B1 = {h1, mb1, h2};

    //""" -- b2 -- """
    dip_nam =  "b2";
    dip_len =  1.228262;
    dip_ang =  4.103510 * deg_2_rad;
    dip_K   = -0.78;
    dip_S   =  0.00;
    h1      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 1*dip_ang/2, 0*dip_ang/2, dip_K, dip_S);
    h2      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 0*dip_ang/2, 1*dip_ang/2, dip_K, dip_S);
    std::vector<Element> B2 = {h1, mb2, h2};

    //""" -- b3 -- """
    dip_nam =  "b3";
    dip_len =  0.428011;
    dip_ang =  1.429950 * deg_2_rad;
    dip_K   = -0.78;
    dip_S   =  0.00;
    h1      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 1*dip_ang/2, 0*dip_ang/2, dip_K, dip_S);
    h2      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 0*dip_ang/2, 1*dip_ang/2, dip_K, dip_S);
    std::vector<Element> B3 = {h1, mb3, h2};

    //""" -- bc -- """
    dip_nam =  "bc";
    dip_len =  0.125394;
    dip_ang =  1.4 * deg_2_rad;
    dip_K   =  0.00;
    dip_S   = -18.93;
    Element bce      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 1*dip_ang/2, 0*dip_ang/2, dip_K, dip_S);
    Element bcs      = Element::rbend(dip_nam, dip_len/2, dip_ang/2, 0*dip_ang/2, 1*dip_ang/2, dip_K, dip_S);
    std::vector<Element> BC = {bce, mc, bcs};

    //""" --- correctors --- """
    Element ch     = Element::hcorrector("hcm",  0, 0);
    Element cv     = Element::vcorrector("vcm",  0, 0);
    Element crhv   = Element::corrector ("crhv", 0, 0, 0);

    //""" --- sextupoles --- """
    Element sa1      = Element::sextupole("sa1", 0.150000, sa1_strength);
    Element sa2      = Element::sextupole("sa2", 0.150000, sa2_strength);
    Element sb1      = Element::sextupole("sb1", 0.150000, sb1_strength);
    Element sb2      = Element::sextupole("sb2", 0.150000, sb2_strength);
    Element sd1      = Element::sextupole("sd1", 0.150000, sd1_strength);
    Element sf1      = Element::sextupole("sf1", 0.150000, sf1_strength);
    Element sd2      = Element::sextupole("sd2", 0.150000, sd2_strength);
    Element sd3      = Element::sextupole("sd3", 0.150000, sd3_strength);
    Element sf2      = Element::sextupole("sf2", 0.150000, sf2_strength);

    //""" --- rf cavity --- """
    Element cav = Element::rfcavity("cav", 0, 500e6, 2.5e6);

    //""" lines """
    std::vector<Element> insa   = { dia1, mida, dia2, crhv, cv, d12, ch, d12, sa2, d12, mon, d12, qaf, d23, qad, d17, sa1, d17};
    std::vector<Element> insb   = { dib1, midb, dib2, d10, crhv, qbd2, d12, cv, d12, ch, d12, sb2, d12, mon, d12, qbf, d23, qbd1, d17, sb1, d17};
    std::vector<Element> cline1 = { d32, cv,  d12, ch,  d15, sd1, d17, qf1, d12, mon, d11, sf1, d20, qf2, d17, sd2, d12, ch, d10, mon, d10};
    std::vector<Element> cline2 = { d18, cv,  d26, sd3, d17, qf3, d12, mon, d11, sf2, d20, qf4, d15, ch,  crhv, d12, mon, d44};
    std::vector<Element> cline3 = { d44, mon, d12, ch,  d15, qf4, d20, sf2, d11, mon, d12, qf3, d17, sd3, d26, cv, crhv, d18};
    std::vector<Element> cline4 = { d20, ch,  d12, sd2, d17, qf2, d20, sf1, d11, mon, d12, qf1, d17, sd1, d15, ch,  d12, cv, d22, mon, d10};

    //""" Injection Section """
    Element dmiainj  = Element::drift("dmiainj", 0.3);
	Element dinjk3   = Element::drift("dinjk3" , 0.3);
	Element dk3k4    = Element::drift("dk3k4"  , 0.6);
	Element dk4pmm   = Element::drift("dk4pmm" , 0.2);
	Element dpmmcv   = Element::drift("dpmmcv" , (3.2692 + 3.65e-3 - 0.3 - 0.3 - 0.6 - 0.2 - 3*0.6));
	Element dcvk1    = Element::drift("dcvk1"  , (3.2692 + 3.65e-3 - 0.6 - 1.4 - 2*0.6));

	Element dk1k2    = Element::drift("dk1k2"  , 0.6);
	Element sef      = Element::sextupole("sef", 0.6, 0.0, 5);
	Element dk2sef   = Element::drift("dk2mia" , 0.8);

	Element kick     = Element::corrector("kick", 0.6, 0, 0);
	Element pmm      = Element::sextupole("pmm", 0.6, 0.0, 5);
	Element inj      = Element::marker("inj");

	std::vector<Element> insaend  = {cv, d12, ch, d12, sa2, d12, mon, d12, qaf, d23, qad, d17, sa1, d17};
    std::vector<Element> insainj  = latt_join({{dmiainj, inj, dinjk3, kick, dk3k4, kick, dk4pmm, pmm, dpmmcv}, insaend});
    std::vector<Element> injinsa  = latt_join({latt_reverse(insaend), {dcvk1, kick, dk1k2, kick, dk2sef, sef}});

    std::vector<Element> B3BCB3   = latt_join({B3,{d13},BC,{d13},B3});
    std::vector<Element> R01 = latt_join({injinsa, {fim, inicio, mia}, insainj});   //#% injection sector, marker of the lattice model starting element
    std::vector<Element> R03 = latt_join({latt_reverse(insa), {mia, cav}, insa});   //#% sector with cavities
    std::vector<Element> R05 = latt_join({latt_reverse(insa), {mia}, insa});
    std::vector<Element> R07(R05);
    std::vector<Element> R09(R05);
    std::vector<Element> R11(R05);
    std::vector<Element> R13(R05);
    std::vector<Element> R15(R05);
    std::vector<Element> R17(R05);
    std::vector<Element> R19(R05);
    std::vector<Element> R02 = latt_join({latt_reverse(insb), {mib}, insb});
    std::vector<Element> R04(R02);
    std::vector<Element> R06(R02);
    std::vector<Element> R08(R02);
    std::vector<Element> R10(R02);
    std::vector<Element> R12(R02);
    std::vector<Element> R14(R02);
    std::vector<Element> R16(R02);
    std::vector<Element> R18(R02);
    std::vector<Element> R20(R02);


    std::vector<Element> C01 = latt_join({B1, cline1, B2, cline2, B3BCB3, cline3, B2, cline4, B1});
    std::vector<Element> C02(C01);
    std::vector<Element> C03(C01);
    std::vector<Element> C04(C01);
    std::vector<Element> C05(C01);
    std::vector<Element> C06(C01);
    std::vector<Element> C07(C01);
    std::vector<Element> C08(C01);
    std::vector<Element> C09(C01);
    std::vector<Element> C10(C01);
    std::vector<Element> C11(C01);
    std::vector<Element> C12(C01);
    std::vector<Element> C13(C01);
    std::vector<Element> C14(C01);
    std::vector<Element> C15(C01);
    std::vector<Element> C16(C01);
    std::vector<Element> C17(C01);
    std::vector<Element> C18(C01);
    std::vector<Element> C19(C01);
    std::vector<Element> C20(C01);


    the_ring = latt_join({
        R01, C01, R02, C02, R03, C03, R04, C04, R05, C05,
        R06, C06, R07, C07, R08, C08, R09, C09, R10, C10,
        R11, C11, R12, C12, R13, C13, R14, C14, R15, C15,
        R16, C16, R17, C17, R18, C18, R19, C19, R20, C20,
    });


    //""" shift lattice to start at the marker "inicio" """
    std::vector<int> idx = latt_findcells_fam_name(the_ring, "inicio");
    if (idx.size() > 0) {
    	std::vector<Element>::iterator it = the_ring.begin() + idx[0];
    	std::rotate(the_ring.begin(), it, the_ring.end());
    };

    //""" check if there are elements with negative lengths """
    std::vector<double> lens = latt_getcellstruct<double>(the_ring, "length", latt_range(the_ring));
    for(unsigned int i=0; i<lens.size(); ++i) {
    	if (lens[i] < 0) {
    		std::cerr << "negative drift in lattice!" << std::endl;
    	}
    }

    //""" sets cavity frequency according to lattice length """
    double C = latt_findspos(the_ring, 1+the_ring.size());

    double rev_freq = light_speed / C;
    std::vector<int> rf_idx = latt_findcells_fam_name(the_ring, "cav");

    for(unsigned int idx = 0; idx<rf_idx.size(); ++idx) {
    	the_ring[rf_idx[idx]].frequency = rev_freq * harmonic_number;
    }
    latt_setcavity(the_ring, "on");
    //latt_setradiation(the_ring, "on", 3e9);

    //""" adjusts number of integraton steps for each element family """
    the_ring = latt_set_num_integ_steps(the_ring);

}
