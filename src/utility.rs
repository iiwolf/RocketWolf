/* misc utility functions */

    pub static DEBUG: bool = true;
    pub static PI:f32 = 3.14159265359;  //pumpkin
    pub static g0_EN: f32 = 32.17405;   //ft/s^2
    pub static g0_SI: f32 = 9.80665;    //m/s^2
    pub static p_atm_EN: f32 = 14.696;  //psi
    pub static p_atm_SI: f32 = 101325.0;  //Pa
    pub static Rprime_EN: f32 = 1544.0; //ft-lbf / lb-mol-R*
    pub static Rprime_SI: f32 = 8314.3; //J/kg-mol-K

    //print result with variable name, value, and units
    pub fn result(variable_name: &str, variable_value: f32, variable_units: &str){
        println!( "{}: {} {}", variable_name, variable_value, variable_units);
    }

    //same as result, but can be turned off
    pub fn debug(variable_name: &str, variable_value: f32, variable_units: &str){
        if DEBUG {
            result(variable_name,variable_value,variable_units);
        } 
    }

    //problem header
    pub fn problem_header(header: &str){
        println!("\n -- Problem {} --", header);

    }

    fn tests(){

    let k = 1.2;
    let R = 360.0;
    let p1 = 2.039e6;
    let T1 = 2800.0;
    let p2_p1 = 1.0/800.0;
    let v2 = f32::sqrt(2.0*k*R*T1 / (k-1.0) 
                * (1.0 - f32::powf(p2_p1,(k-1.0)/k)));
    debug("v2",v2,"ms/s"); 
    }
