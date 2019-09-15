use crate::utility;

/*
    FIND:
        v2, T2, CF for an optimum exapnsion nozzle at sea level
    Equations:
        * Perfect Gas Law: pV = RT
        * R = R' / M
*/
fn sutton3_2(){ 

    utility::problem_header("3-2");

    //GIVEN
    let k = 1.3;        //-
    let mdot = 3.7;     //kg/s
    let p1 = 2.1e6;     //Pa
    let T1 = 2585.0;    //K
    let M = 18.0;       //kg/kg-mol

    //since we are at sea level + optimum expansion
    let p2 = utility::p_atm_SI;

    //find R for given M
    let R = utility::Rprime_SI / M;
    utility::debug("R",R,"");

    //v2
    let v2 = f32::sqrt(2.0*k*R*T1 / (k-1.0) 
                * (1.0 - f32::powf(p2/p1,(k-1.0)/k)));
    utility::result("v2",v2,"m/s");

    //get T2 via isentropic relationships
    let T2 = T1 * f32::powf(p2/p1,1.0 - 1.0/k);
    utility::result("T2",T2,"K");

    //get area ratio from 3-4 or eq 3-25
    let At_A2 = f32::powf((k + 1.0)/2.0,1.0/(k-1.0)) * f32::powf(p2/p1,1.0/k)
                * f32::sqrt((k + 1.0) / (k - 1.0) * (1.0 - f32::powf(p2/p1,(k-1.0)/k)));
    utility::debug("At_A2",At_A2,"");

    //using ideal gas equation with mass flow to get At and subsequently A2
    let At = mdot/p1 * f32::sqrt(R*T1 / (k * f32::powf(2.0 / (k + 1.0),(k + 1.0)/(k - 1.0))));
    let A2 = At / At_A2;
    utility::debug("A2",A2,"m");    

    //thrust coefficient with p3 = p2 so last two terms cancel
    let V2 = R*T2 / p2;     //first get V2 from ideal gas equation
    let CF = f32::powf(v2,2.0) * A2 / (p1 * At * V2);
    utility::result("CF",CF,"");

    //  -- Problem 3-2 --
    // R: 461.90555
    // v2: 2281.9084 m/s
    // T2: 1284.2579 K
    // At_A2: 0.30394876 
    // A2: 0.009492756 m
    // CF: 1.3934388 
}

/*
    FIND:
        (a) vt
        (b) Vt
        (c) mdot and Is
        (d) Thrust (F)
        (e) Mach at throat

    EQUATIONS:
        Specific volume at the throat:
            Vt = V1*f32::pow((k + 1.0) / 2.0, 1 / (k - 1))
        Temperature at throat
            Tt = 2.0 * T1 / (k + 1.0)
        Velocity at throat:
            vt = sqrt(k * R * Tt) = at
        Mass flow from conservation of mass
            mdot = A v / V
    
    ASSUMPTIONS:
        Ideal rocket
        Isentropic + no heat losses

    COMMENTS:
        _Ideal_ implies no losses, whereas _optimum_ is a separate concept reflecting he best
        calculated performace at a particular set of given pressures.
*/
fn sutton3_5() -> (f32, f32, f32, f32) {

    utility::problem_header("3-5");
    
    //GIVEN
    let M = 24.0;           //kg/kg-mol
    let p1 = 2.533e6;       //Pa
    let p2 = 0.090e6;       //Pa
    let T1 = 2900.0;        //K
    let At = 0.00050;       //m^2
    let k = 1.3;            //-

    //preliminary calcs
    let R = utility::Rprime_SI / M;

    //get Tt to plug into vt
    let Tt =  2.0 * T1 / (k + 1.0);
    let vt = f32::sqrt(k * R * Tt);
    utility::result("vt", vt, "m/s");

    //get V1 from ideal gas equations to get Vt
    let V1 = R*T1 / p1; 
    let Vt = V1*f32::powf((k + 1.0) / 2.0, 1.0 / (k - 1.0));
    utility::result("Vt", Vt, "m^3 / kg");

    //get mdot from conservation of mass
    let mdot = At * vt / Vt;
    utility::result("mdot", mdot, "kg / s");

    //thrust from ideal rocket equation with constant k and p2 = p3 (3-29)
    let F = At*p1*f32::sqrt(
        2.0*f32::powf(k,2.0)/(k - 1.0) 
        * f32::powf(2.0/(k + 1.0),(k+1.0)/(k-1.0)) 
        * (1.0 - f32::powf(p2/p1,(k - 1.0)/k)));
    utility::result("F", F, "N");

    //with thrust and mdot, get specific impulse
    let Is = F / (mdot*utility::g0_SI);
    utility::result("Is", Is, "s");

    //mach number at throat 
    let Mt = vt / f32::sqrt(k*R*Tt);
    utility::result("Mt", Mt, "");

    //  -- Problem 3-5 --
    // mdot: 0.8431321 kg / s
    // F: 1823.2208 N
    // Is: 220.50731 s
    // Mt: 1

    //return values for problem 6
    return (p1, p2, At, F);
}


/*
    GIVEN:
        Results from problem 5
    FIND:
        Determine the ideal thrust from problem 5 by two methods    
    ASSUME:
        Ideal rocket
        Isentropic
*/
fn sutton3_6(p1: f32, p2: f32, At: f32, F:f32){

    utility::problem_header("3-6");

    //via thrust coefficient equation 3-30
    let CF_equation = F/(p1 * At);
    utility::result("CF (equation)", CF_equation, "");

    //get pressure ratio to use in figure 3-5 table lookup
    let p1_p2 = p1 / p2;
    utility::result("p1/p2", p1_p2, "");

    let CF_table = 1.45;
    utility::result("CF (table)", CF_table, "");
    
    //  -- Problem 3-6 --
    // CF (equation): 1.4395742
    // p1/p2: 28.144444
    // CF (table): 1.45

}


/*
    GIVEN:
        Ideal rocket with given characteristics
    FIND:
        (a) chamber pressure (p1)
        (b) thrust (F)
        (c) specific impulse (Is)
    ASSUMPTIONS:
        Ideal Rocket
        Isentropic
    EQUATIONS:
        c_star = p1 * At / mdot = Is * g0 / CF
        F = CF * At * p1
        Is = F / mdot * g0
*/
fn sutton3_8(){

    utility::problem_header("3-8");

    //given
    let c_star = 1500.0;    //m/s
    let At = 20.0 / 100.0;  //m
    let CF = 1.38;          //-
    let mdot = 40.0;        //kg/s

    //use characterstic velocity eq to get p1
    let p1 = c_star * mdot / At;
    utility::result("p1",p1,"Pa");

    //with p1, we are able to get force
    let F = CF * At * p1;
    utility::result("F",F,"N");

    //specific impulse by characteristic velocity or thrust eq
    let Is = c_star * CF / utility::g0_SI;
    utility::result("Is",Is,"s");
    
    let Is = F / (mdot * utility::g0_SI);
    utility::result("Is",Is,"s");

    //  -- Problem 3-8 --
    // p1: 300000 Pa
    // F: 82800 N
    // Is: 211.08125 s
    // Is: 211.08127 s
}

/*
    GIVEN:
        Design a supersonic nozzle at the given conditions
    FIND:
        (a) M2, v2, T2
        (b) p1
        (c) What happens to thrust and exit velocity if p1 is doubled
        (d) How close to optimum nozzle expansion is this nozzle?
    ASSUME:
        Ideal rocket
        Isentropic
        No change in gas properties for part 2
*/
fn sutton3_12(){

    utility::problem_header("3-12");

    //GIVEN
    let alt = 10.0e3;       //m
    let A2_At = 8.0;        //-
    let T0 = 3000.0;        //K
    let R = 378.0;          //J/kg-K
    let k = 1.3;            //-

    //p_atm at 10km
    let p3 = 26436.27;      //Pa
    let p2 = p3;            //for optimum conditions

    //isentropic relationships given (k,A/At) for M2, T2/T0, and P2/P0
    let M2 = 3.385;
    let T2_T0 = 0.3677;
    let p2_p0 = 0.0131;

    //use ratios and givens to solve for unknowns
    let T2 = T2_T0 * T0;
    let p0 = p2 / p2_p0;

    //v2 from definition of mach number
    let v2 = M2 * f32::sqrt(k*R*T2);

    utility::result("M2", M2, "");
    utility::result("v2", v2, "m/s");
    utility::result("T2", T2, "K");

    //chamber properties are equal to stagnation
    let p1 = p0;
    let T1 = T0;
    utility::result("p1", p1, "Pa");

    //v2 sanity check to see if T1 and p1 are really correct
    let v2 = f32::sqrt(2.0*k*R*T1 / (k-1.0) 
                * (1.0 - f32::powf(p2/p1,(k-1.0)/k)));
    utility::debug("v2",v2,"m/s");

    //  -- Problem 3-12 --
    // M2: 3.385
    // v2: 2492.2046 m/s
    // T2: 1103.1 K
    // p1: 2018035.9 Pa
    // v2: 2492.7878 m/s
}

pub fn main() {
    println!("\n ~~ MAE540 - HW03 ~~ \n");  
    sutton3_2();
    let (p1, p2, At, F) = sutton3_5();
    sutton3_6(p1,p2,At,F);
    sutton3_8();
    sutton3_12();
}




