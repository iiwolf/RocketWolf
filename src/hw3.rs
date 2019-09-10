use crate::utility;


/*
Qs:
    Is temperature isentropic across throat/shock?? No right? How is the T2/T1 off of
    P2/P1 and k arrived at

*/
pub fn main() {
    println!("\n ~~ MAE540 - HW03 ~~ \n");  
    sutton3_2();
    let (p1, p2, At, F) = sutton3_5();
    sutton3_6(p1,p2,At,F);
    sutton3_8();
    sutton3_12();
}


/*
    Find v2, T2, CF for an optimum exapnsion nozzle at sea level

    Equations:
        * Perfect Gas Law: pV = RT
        * R = R' / M
*/
fn sutton3_2(){ 

    utility::problem_header("3-2");

    //given
    let k = 1.3;        //-
    let mdot = 3.7;     //kg/s
    let p1 = 2.1e6;     //Pa
    let T1 = 2585.0;      //K
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

}

/*
    Given the optimum rocket propulsion system
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
        No heat losses

    COMMENTS:
        _Ideal_ implies no losses, whereas _optimum_ is a separate concept reflecting he best
        calculated performace at a particular set of given pressures.
*/
fn sutton3_5() -> (f32, f32, f32, f32) {

    utility::problem_header("3-5");
    
    //given
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

    //thrust from ideal rocket equation with constant k and p2 = p3
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

    //return values for problem 6
    return (p1, p2, At, F);
}


/*
    Determine the ideal thrust from problem 5 by two methods
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
        No change in gas properties
        Isentropic
*/
fn sutton3_12(){

    utility::problem_header("3-12");

    //given
    let alt = 10.0e3;       //m
    let A2_At = 8.0;        //-
    let T0 = 3000.0;        //K
    let R = 378.0;          //J/kg-K
    let k = 1.3;            //-

    //p_atm at 10km
    let p3 = 26436.27;      //Pa
    let p2 = p3;            //for optimum conditions
    let p2_p1 = f32::powf(1.0 + (k - 1.0)/2.0*f32::powf(M2,2.0),-k/(K-1.0));
    //isentropic relationships given (k,A/At) for M2, T2/T0, and P2/P0
    let M2 = 3.385;
    let T2_T0 = 0.3677;
    // let T2_Tt = 0.4229;
    let p2_p0 = 0.0131;
    // let p2_pt = 0.0240;

    //use ratios and givens to solve for unknowns
    let T2 = T2_T0 * T0;
    let p0 = p2 / p2_p0;
    // let Tt = T2 / T2_Tt;

    //v2 from definition of mach number
    let v2 = M2 * f32::sqrt(k*R*T2);

    utility::result("M2", M2, "");
    utility::result("v2", v2, "m/s");
    utility::result("T2", T2, "K");

    //chamber properties are equal to stagnation
    let p1 = p0;
    let T1 = T0;

    //get pressure at throat to use in eq 3-20 
    // let pt = p2 / p2_pt;
    // let p1 = pt / f32::powf(2.0 / (k + 1.0), k / (k - 1.0));
    utility::result("p1", p1, "Pa");

    //v2 sanity check to see if T1 and p1 are really correct
    let v2 = f32::sqrt(2.0*k*R*T1 / (k-1.0) 
                * (1.0 - f32::powf(p2/p1,(k-1.0)/k)));
    utility::debug("v2",v2,"m/s");
    //get T1 to check above result with v2 equation
    // let T1 = Tt / (2.0*(k + 1.0));
    // utility::debug("T1",T1,"K");

    //what happens to thrust + velocity if pressure is doubled (but no other changes in gas properties)
    let p1 = p1 * 2.0;
    let v2 = f32::sqrt(2.0*k*R*T1 / (k-1.0) 
                * (1.0 - f32::powf(p2/p1,(k-1.0)/k)));

    //thrust from ideal rocket equation with constant k and p2 = p3
    let F = At*p1*f32::sqrt(
        2.0*f32::powf(k,2.0)/(k - 1.0) 
        * f32::powf(2.0/(k + 1.0),(k+1.0)/(k-1.0)) 
        * (1.0 - f32::powf(p2/p1,(k - 1.0)/k)));

    utility::result("F", F, "N");
    utility::result("v2",v2,"m/s"); 
}

fn A2_At(M2: f32, k: f32) -> (f32){
    return (1.0/M2) * f32::powf(
        (2.0/(k+1.0)) * (1.0 + (k - 1.0)/2.0*f32::powf(M2,2.0)),(k + 1.0)/(2.0*(k-1.0)));
}

// fn error_mach(A2_At: f32, M2i: f32, k: f32){
//     return 
// }
fn SP03(){


}