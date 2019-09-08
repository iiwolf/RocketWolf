use gnuplot::*;

static DEBUG: bool = false;
static g0_EN: f32 = 32.17405;   //ft/s^2
static g0_SI: f32 = 9.80665;    //m/s^2
static P_atm: f32 = 14.696;     //psi
static PI:f32 = 3.14159265359;  //pumpkin

fn main() {
    println!("\n ~~ MAE540 - HW02 ~~ \n");
    sutton2_3();
    sutton2_7();
    sutton2_8();
    SP02B();
}

/*
Sutton 2-3
A certain rocket engine (flying horizontally) has an effective exhaust velocity of
7000 ft/sec; it consumes 280 lbm/sec of propellant mass, and liberates 2400 Btu/lbm.
The unit operates for 65 sec. Construct a set of curves plotting the propulsive, internal,
and overall efficiencies versus the velocity ratio u∕c (0 < u∕c < 1.0). The rated flight
velocity equals 5000 ft/sec. Calculate 
    (a) the specific impulse;
    (b) the total impulse;
    (c) the mass of propellants required;
    (d) the volume that the propellants occupy if their average specific gravity is 0.925. Neglect gravity and drag.
Answer: (a) 217.4 sec; (b) 3,960,000 lbf-sec; (c) 18,200 lbm; (d) 315 ft3.
*/

fn sutton2_3(){

    problem_header("2-3");

    //given
    let c = 7000.0;     //ft/s
    let mdot = 280.0;   //lbm/s
    let Qr = 2400.0;    //BTU/lbm
    let t = 65.0;       //s
    let v = 5000.0;     //ft/s
    let SG = 0.925;     //unitless - lbm / ft^3

    //constants
    let g0 = g0_EN;         //ft/s^2
    let rho_water = 62.4;   //lbs/ft
    let J = 778.0;          //ft-lbf/BTU

    //thrust from mass flow and exhaust velocity
    let F = mdot * c;
    let Is = specific_impulse(F,mdot,g0);  //seconds
    result("Specific Impulse",Is,"s"); 

    //total impulse given specific impulse, mass propellant
    let It = lbm_to_lbf(F*t);
    result("Total Impulse",It,"lbf-s");

    //total mass of propellant required
    let mtot = mdot * t;
    result("Mass Propellant", mtot, "lbm");

    //total volume occupied by propellant
    let rho_propellant = SG*rho_water;
    let V = mtot / rho_propellant;
    result("Volume of Propellant", V, "ft^3");

    let p_chem = mdot * Qr * J;
    let n_comb = 0.95;
    let mut n_int = Vec::new(); //interal efficiency
    let mut n_p = Vec::new();   //propulsive efficiency
    let mut uc = Vec::new();    //uc step
    let num_it = c as i32;
    let step_size = 1000;

    //loop over velocities, stepping by step_size
    for u in (0..=num_it).step_by(step_size) {
        uc.push(u as f32 / c);
        n_int.push(internal_efficiency(mdot / g0_EN, u as f32, n_comb, p_chem));
        n_p.push(propulsive_efficiency(u as f32,c));
    }

    //plot u/c from 0 to 1
    let mut fg = Figure::new();
    fg.axes2d()
    .lines(
        &uc,
		&n_int,
		&[Caption("Internal Efficiency")],
    );

    let mut fg2 = Figure::new();    
    fg.axes2d()
    .lines(
        &uc,
		&n_p,
		&[Caption("Propulsive Efficiency")],
    );
    
    fg.show();    

}

/* 
Sutton 2-7
For a solid propellant rocket motor with a sea-level thrust of 207,000 lbf, determine:
FIND:
    (a) the (constant) propellant mass flow rate ṁ and the specific impulse Is at sea level,
    (b) the altitude for optimum nozzle expansion as well as the thrust and specific impulse 
        at this optimum condition and
    (c) at vacuum conditions.

GIVEN:
The initial total mass of the rocket motor is 50,000 lbm and its propellant mass fraction is 0.90. 
The residual propellant (called slivers, combustion stops when the chamber pressure falls below a deflagration limit) amounts
to 3% of the burnt. The burn time is 50 seconds; the nozzle throat area (At) is 164.2 in.2
and its area ratio (A2/At) is 10. The chamber pressure (p1) is 780 psia and the pressure ratio
(p1/p2) across the nozzle may be taken as 90.0. Neglect any start/stop transients and use
the information in Appendix 2.

Answers: (a) ṁ = 873 lbm/sec, 237 sec., (b) F = 216,900 lbf, Is = 248.5 sec., (c) F =
231,000 lbf, Is = 295 sec.
*/

fn sutton2_7(){

    problem_header("2-7");

    //given
    let rocket_thrust = 207000.0;   //lbf
    let rocket_mass = 50000.0;      //lbm
    let prop_mass_fraction = 0.90;  //- 
    let residual_propellant = 0.03; //%
    let t_burn = 50.0;              //s
    let At = 164.2;                 //in^2
    let A2_At = 10.0;               //-
    let P1 = 780.0;                 //psia
    let P1_P2 = 90.0;               //-

    //(a)
    let propellant_mass = rocket_mass * prop_mass_fraction;
    let mdot = (1.0 - residual_propellant) * propellant_mass / t_burn;
    result("Propellant Mass Flow",mdot,"lbm/s");

    let Is = specific_impulse(rocket_thrust,mdot,1.0);
    result("Specific Impulse", Is, "s");

    //(b)
    let A2 = A2_At * At;
    let P2 = P1 / P1_P2;
    let pressure_ratio = P2 / P_atm;

    //look up altitude based on ratio
    let alt = 5e3;
    result("Altitude",alt,"ft");

    let v2 = (rocket_thrust - (P2 - P_atm) * A2) / mdot;
    let F_optimum = mdot * v2;
    result("F_optimum",F_optimum,"lbf");
     
    let Is = specific_impulse(F_optimum,mdot,1.0);
    result("Specific Impulse",Is,"s");

    let F_vacuum = mdot * v2 + P2 * A2;
    result("F_vacuum", F_vacuum, "lbf");

    let Is = specific_impulse(F_vacuum,mdot,1.0);
    result("Specific Impulse",Is,"s");       
}

/*
Sutton 2-8
During the boost phase of the Atlas V, the RD-180 engine operates together with three
solid propellant rocket motors (SRBs) for the initial stage. For the remaining thrust time,
the RD-180 operates alone. Using the information given in Table 1–3, calculate the overall
effective exhaust velocity for the vehicle during the initial combined thrust operation
Answer: 309s (?) 9,495 ft/s
*/

fn sutton2_8(){

    problem_header("2-8");

    //specific impulses of booster and RD-180 engine given by table 1-3
    let Is_solid_booster = 279.3;   //(s) at sea level
    let Is_RD_180 = 310.7;          //(s) at sea level
    let g0 = g0_EN;                 //(ft/s^2)

    //problem states that Atlas works with 3 boosters and 1 engine for 
    // the initial combined thrust operation
    let c = (3.0*Is_solid_booster + 1.0*Is_RD_180)*g0;

    result("Effective Exhaust Velocity", c, "ft/s");
}

/*
SP02B
An Ideal Rocket motor operating with a converging-diverging nozzle. The chamber
pressure, P1 is 12 MPa, the chamber temperature, T1 is 3,000K, the specific heat ratio, k, of the
gases is 1.3. The cross section of the throat and the exit flow diameter are shown in the schematic
below. Clearly define the control volume and the coordinate system used for the analysis. The
ambient pressure, P3, is 1 atmosphere. Specify throat properties with subscript t.
a. The thrust of the rocket motor (set up equations only).
b. The axial stress at the throat (set up equations only).
c. Does the axial stress represent tension or compression?
d. Bonus: Calculate the axial stress
*/

fn SP02B(){
    
    problem_header("SP02B");

    //given
    let P1 = 12.0e6;        //Pa
    let P3 = 101325.0;      //Pa
    let T1 = 3000;          //K
    let k = 1.3;            //
    let rt = 0.300/2.0;     //m
    let r2 = 0.600/2.0;     //m
    let thickness = 0.006;  //m

    //get areas
    let At = PI*(f32::powf(rt,2.0) - f32::powf(rt - thickness,2.0));
    let A2 = PI*(f32::powf(r2,2.0));
    let A2_At = A2 / At;
    result("A2/At", A2_At, " - ");

    //then from figure 3-5 we get
    let P1_P2 = 900.0;
    let P2 = P1/P1_P2;

    //which is the final piece we need for the ideal thrust equation (3-29)
    let F = At*P1*f32::sqrt(
            2.0*f32::powf(k,2.0)/(k - 1.0) 
            * f32::powf(2.0/(k + 1.0),(k+1.0)/(k-1.0)) 
            * (1.0 - f32::powf(P2/P1,(k - 1.0)/k)))
            + (P2 - P3)*A2;
    result("Thrust", F, "N");

    //which allows us to calculate stress on area At
    let sigma = F / At;
    result("Stress", sigma, "N / m^2");

}

//calculates impulse for single time step
fn impulse(thrust_force: f32, timestep: f32) -> f32
{
    return thrust_force*timestep;
}

//specific impulse: given total impulse, mass of propellant, and g0
fn specific_impulse(total_impulse: f32, mass_propellant_expelled: f32, acceleration_of_gravity: f32) -> f32{
    return total_impulse / (mass_propellant_expelled * acceleration_of_gravity);
}

//effective exhaust velocity given specific impulse and g0
fn effective_exhaust_velocity(specific_impulse: f32, acceleration_of_gravity: f32) -> f32{
    return specific_impulse*acceleration_of_gravity;
}

//convert pounds mass to pounds force
fn lbm_to_lbf(value_in_lbm: f32) -> f32{
    return value_in_lbm / g0_EN;
}

//internal efficiency given mass flow, velocity, combustive efficiency and pchem
fn internal_efficiency(mdot: f32, u: f32, n_comb: f32, p_chem: f32) -> f32{
    return 0.5 * mdot * f32::powi(u,2) / (n_comb * p_chem);
}

//propulsive efficient given u and c
fn propulsive_efficiency(u: f32, c: f32) -> f32{
    return 2.0*(u/c) / (1.0 + f32::powi(u/c,2));
}

//thrust
fn thrust(mdot: f32, v2: f32, p2: f32, p3: f32, A2: f32) -> f32{
    return mdot * v2 + (p2 - p3) * A2;
}

/* UTILITY FUNCTIONS */

fn result(variable_name: &str, variable_value: f32, variable_units: &str){
    println!( "{}: {} {}", variable_name, variable_value, variable_units);
}

fn debug(variable_name: &str, variable_value: f32, variable_units: &str){
    if DEBUG {
        result(variable_name,variable_value,variable_units);
    } 
}

//problem header
fn problem_header(header: &str){
    println!("\n -- Problem {} --", header);
}