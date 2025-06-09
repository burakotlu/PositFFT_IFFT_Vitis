#include "posit.hpp"
#include <cmath>
#include <hls_math.h>


double d_Angle_Array[MAX_SIZE] = {0.0};
size_t d_Angle_Counter = 0;

void addTod_Angle_Array(double value) {
    if (d_Angle_Counter < MAX_SIZE) {
        d_Angle_Array[d_Angle_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
float f_Angle_Array[MAX_SIZE] = {0.0};
size_t f_Angle_Counter = 0;

void f_addTof_Angle_Array(float value) {
    if (f_Angle_Counter < MAX_SIZE) {
        f_Angle_Array[f_Angle_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
double p_Angle_Array[MAX_SIZE] = {0.0};
size_t p_Angle_Counter = 0;

void p_addTop_Angle_Array(double value) {
    if (p_Angle_Counter < MAX_SIZE) {
        p_Angle_Array[p_Angle_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
////////////////////////////////////////////
double d_RP_Array[MAX_SIZE] = {0.0};
size_t d_RP_Counter = 0;

void addTod_RP_Array(double value) {
    if (d_RP_Counter < MAX_SIZE) {
        d_RP_Array[d_RP_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
float f_RP_Array[MAX_SIZE] = {0.0};
size_t f_RP_Counter = 0;

void f_addTof_RP_Array(float value) {
    if (f_RP_Counter < MAX_SIZE) {
        f_RP_Array[f_RP_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
double p_RP_Array[MAX_SIZE] = {0.0};
size_t p_RP_Counter = 0;

void p_addTop_RP_Array(double value) {
    if (p_RP_Counter < MAX_SIZE) {
        p_RP_Array[p_RP_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
//////////
double d_IMG_Array[MAX_SIZE] = {0.0};
size_t d_IMG_Counter = 0;

void addTod_IMG_Array(double value) {
    if (d_IMG_Counter < MAX_SIZE) {
        d_IMG_Array[d_IMG_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
float f_IMG_Array[MAX_SIZE] = {0.0};
size_t f_IMG_Counter = 0;

void f_addTof_IMG_Array(float value) {
    if (f_IMG_Counter < MAX_SIZE) {
        f_IMG_Array[f_IMG_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
double p_IMG_Array[MAX_SIZE] = {0.0};
size_t p_IMG_Counter = 0;

void p_addTop_IMG_Array(double value) {
    if (p_IMG_Counter < MAX_SIZE) {
        p_IMG_Array[p_IMG_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
//////////
double d_DTH_Array[DTHETA_SIZE] = {0.0};
size_t d_DTH_Counter = 0;

void addTod_DTH_Array(double value) {
    if (d_DTH_Counter < DTHETA_SIZE) {
        d_DTH_Array[d_DTH_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
float f_DTH_Array[DTHETA_SIZE] = {0.0};
size_t f_DTH_Counter = 0;

void f_addTof_DTH_Array(float value) {
    if (f_DTH_Counter < DTHETA_SIZE) {
        f_DTH_Array[f_DTH_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}
double p_DTH_Array[DTHETA_SIZE] = {0.0};
size_t p_DTH_Counter = 0;

void p_addTop_DTH_Array(double value) {
    if (p_DTH_Counter < DTHETA_SIZE) {
        p_DTH_Array[p_DTH_Counter++] = value;
    } else {
        std::cerr << "Global array overflow!" << std::endl;
    }
}

bool isGreater_posit(ps_t x, ps_t y) {
    // Extract fields
    m_add_t x_mantissa = x.mantissa;
    m_add_t y_mantissa = y.mantissa;

    bool x_sign = x.sign;
    bool y_sign = y.sign;

    regime_t x_regime = x.regime;
    regime_t y_regime = y.regime;

    exponent_t x_exponent = x.exponent;
    exponent_t y_exponent = y.exponent;

    if (x_sign != y_sign)
        return (y_sign == 1);

    bool absXGreaterEqual = 
        (x_regime > y_regime) ||
        (x_regime == y_regime && x_exponent > y_exponent) ||
        (x_regime == y_regime && x_exponent == y_exponent && x_mantissa > y_mantissa);

    if (x_sign == 0)
        return absXGreaterEqual; 
    else
        return !absXGreaterEqual; 
}

ps_t positMod(ps_t x, ps_t y){
    ps_t quotient,c_x, inter,res;
    sf_t sf=0;
   
    bool small;
    bool x_sign = x.sign;
    bool y_sign = y.sign;
    c_x =x;
    c_x.sign = false;

    small= isGreater_posit(y, c_x);
    if (small) return x;
    else{
        quotient = positDiv(c_x,y) ;
        sf = ((sf_t)quotient.regime<<ES)+quotient.exponent;
        quotient.mantissa.range(FRAC_LEN-sf-2,0)=0;
        inter= positMul(quotient,y);
        res = positSub(c_x,inter);
        res.sign = x.sign;

        return res;
    }

}

float floatMul(float x, float y){
    return x*y;    
}

float floatDiv(float x, float y){
    return x/y;    
}
float floatAdd(float x, float y){
    return x+y;
}
float floatSub(float x, float y){
    return x-y;
}
double doubleMul(double x, double y){
    return x*y;
}
double doubleDiv(double x, double y){
    return x/y;
}
double doubleAdd(double x, double y){
    return x+y;
}
double doubleSub(double x, double y){
    return x-y;
}
regime_t LOD(reg_t reg){
	regime_t regime=0;
	bool flag=false;
	bool start=reg.bit(N-2);
	if(start==0) reg=~reg;

	LOD_LOOP:for(int i=N-3;i>=0;i--){
		if(reg.bit(i) && flag==0)
			regime++;
		else flag=true;

	}
    //std::cout<<"regime:"<<regime<<std::endl;
	return start?regime:(~regime)++;
}
sf_t LOD_INT(log_sf_t in){
    bool flag=0;
    int count =1;
    for(int i=in.width-1;i>=0;i-- ){
        if(in[i] ==0){
            if(flag==false) count+=1;
        }
        else flag=true;
    }
    return LOG_IN_SIZE_SF-count;
}
int LOD_ADD(m_add_t in){
    bool flag=0;
    int count =0;
    for(int i=in.width-1;i>=0;i-- ){
        if(in[i] ==0){
            if(flag==false) count+=1;
        }
        else flag=true;
    }
    return count-1;
}
int LOD_MUL(mul_t in){
    bool flag=0;
    int count =0;
    for(int i=in.width-1;i>=0;i-- ){
        if(in[i] ==0){
            if(flag==false) count+=1;
        }
        else flag=true;
    }
    return FRAC_LEN-1-count;
}
int LOD_DIV(dv_t in){
    bool flag=0;
    int count =0;
    for(int i=in.width-1;i>=0;i-- ){
        if(in[i] ==0){
            if(flag==false) count+=1;
        }
        else flag=true;
    }
    return FRAC_LEN-1-count;
}
int LOD_FFT(mantissa_t in){
    bool flag=0;
    int count =0;
    for(int i=in.width-1;i>=0;i-- ){
        if(in[i] ==0){
            if(flag==false) count+=1;
        }
        else flag=true;
    }
    return count;
}
ps_t  decode(posit_t posit){

	ps_t pos;
	bool sign=0, isZero=0, isInf=0;
	regime_t regime;
	exponent_t exponent=0;
	mantissa_t mantissa=0;
	reg_t SREG,REM,reg;

	//decode sign
	sign=posit.bit(N-1);
    //zero case
	if(posit==0)  return pos;
    //inf case
    if(posit == (1<<(N-1))){
        pos.isZero=0;
        pos.isInf = 1;
        return pos;
    }
    //other cases

    if(sign)  posit=~posit+1;
    //find regime
    regime=LOD(posit.range(N-2,0));
    //std::cout<<"regime:"<<regime<<std::endl;
    //find SREG & REM
    if (regime >= 0) {
        if (regime + 3 < N) {
            SREG = (reg_t)regime + 3;
        } else {
            SREG = (reg_t)N;
        }
    } else {            
        SREG = 2 - regime;
    }
    REM=N-SREG;
    //find exponent & mantissa
    if(REM>0){
        if(ES==0){
            mantissa=posit.range(REM-1,0);
        }
        else{
            if(REM>ES){
                exponent=posit.range(REM-1,REM-ES);
                mantissa=posit.range(REM-(ES+1),0);
            }
            else {mantissa=1<<(FRAC_LEN-1); exponent=posit.range(REM-1,0);}
        }
    }
	

	pos.sign=sign;
	pos.regime=regime;
	pos.exponent=exponent;
	pos.mantissa=mantissa;
	pos.isZero=isZero;
    pos.isInf = isInf;
	return pos;
}

posit_t encode(ps_t x){
	posit_t posit=0;
	reg_t SREG,REM,reg;
	regime_t regime;
	exponent_t exponent;
	mantissa_t mantissa;
	bool sign, isZero, isInf;

	isZero=x.isZero;
    isInf = x.isInf;
	sign=x.sign;
	exponent=x.exponent;
	regime=x.regime;
	mantissa=x.mantissa;
    
	if(isZero) posit=0;
    else if(isInf) posit =  1<<(N-1);
    else{
        if (regime >= 0) {
            if (regime + 3 < N) {
                SREG = (reg_t)regime + 3;
            } else {
                SREG = (reg_t)N;
            }
        } else {            
            SREG = 2 - regime;
        }
        REM=N-SREG;
        //SET REGIME BITS
        //max regime case
        if(regime==N-2){
            posit=(~posit).range(N-2,0);
        }
        //positive regime case
        else if(regime>=0){ 
            posit.range(N-2,REM+1)=~posit.range(N-2,REM+1);
        }
        //negative regime case
        else{
            posit.bit(REM)=1;
        }
        //SET OTHER BITS
        if(REM>0){
            if(ES==0){
                posit.range(REM-1,0)=mantissa;
            }
            else{
                if(REM>ES){
                    posit.range(REM-1,REM-ES)=exponent;
                    posit.range(REM-(ES+1),0)=mantissa;
                }
                else posit.range(REM-1,0)=exponent;
            }
        }
	    if(sign) posit=~posit+1;
    }

	return posit;
}

ps_t posit_negate(ps_t posit) {
    ps_t negated_posit;

    // Negate the sign
    negated_posit.sign = (posit.sign == 0) ? 1 : 0;

    // Copy other fields
    negated_posit.regime = posit.regime;
    negated_posit.exponent = posit.exponent; // Assuming "e" refers to exponent
    negated_posit.mantissa = posit.mantissa; // Assuming "m" refers to mantissa
    negated_posit.isZero = posit.isZero;
    negated_posit.isInf = posit.isInf;


    return negated_posit;
}
ps_t float2posit(float in) {
    ps_t result;
    regime_t regime = 0;
    exponent_t exponent = 0;
    mantissa_t mantissa = 0;
    bool sign = false, isZero = false;
    reg_t exact, mant_part, SREG, REM;
    sf_t sf;

    double MAX = pow(2, REG_SHIFT * (N - 2)); //(1 << (REG_SHIFT * (N - 2)));

    // Handle zero input
    if (in == 0.0) {
        isZero = true;
    }

    else if (in >= MAX) {
        regime = N - 2; // Maximum positive regime
    } else if (in <= -MAX) {
        regime = N - 2; // Maximum negative regime
        sign = 1;
    } else {
        // Handle negative numbers
        if (in < 0) {
            sign = true;
            in = -in;
        }


        sf = (int)log2(in); // Scale factor
        exact = 1 << sf;    // Closest power of 2 less than or equal to `in`
        mant_part = (in - exact) * (1 << FRAC_LEN); // Remaining part scaled to mantissa bits
        regime = sf >> ES;

        //find SREG & REM
        if (regime >= 0) {
            if (regime + 3 < N) {
                SREG = (reg_t)regime + 3;
            } else {
                SREG = (reg_t)N;
            }
        } else {            
            SREG = 2 - regime;
        }
        REM=N-SREG;

        if (ES == 0) {
            if (REM > 0) {
                mantissa.range(FRAC_LEN - 2, FRAC_LEN - REM - 1) = mant_part << (REM - sf);
                mantissa.set(FRAC_LEN - 1);
            }
        } else {
            if (REM > ES) {
                exponent = sf & ((1 << ES) - 1); // Extract `ES` bits from sf
                mantissa.range(FRAC_LEN - 2, FRAC_LEN - REM + ES - 1) = mant_part << (REM - sf - ES);
                mantissa.set(FRAC_LEN - 1);
            } else if (REM > 0) {
                exponent = sf & ((1 << REM) - 1); // Use remaining bits for exponent
            }
        }
    }

    result.sign = sign;
    result.regime = regime;
    result.exponent = exponent;
    result.mantissa = mantissa;
    result.isZero = isZero;

    return result;
}
double stable_floor(double x, double epsilon = 1e-9) {
    if (std::abs(x - std::round(x)) < epsilon) {
        return std::round(x);  // Close enough to be considered an integer
    }
    return std::floor(x);
}
ps_t double2posit(double in) {
    ps_t result;
    regime_t regime = 0;
    exponent_t exponent = 0;
    mantissa_t mantissa = 0;

    bool sign = false, isZero = false;
    double exact=0,diff=0;
    reg_t   SREG=0, REM=0;
    mantissa_t mant_part=0;
    mantissa_sf_t mant_with_sf=0, factor=0;
    sf_t sf=0;
    if (in==0){
        result.isZero =1;
        return result;
    }

    if (in < 0) {
        sign = true;
        in = -in;
    }
    
    double sf_d = hls::log2(in);
    double fl= stable_floor(sf_d);
    bool sf_exact =  (fl  == sf_d);  
         
    sf = (sf_t)static_cast<int>(fl); // Scale factor
     
    exact = hls::pow((double)2.0,(double)sf);  
      
    diff = in-exact;
    
    factor.set(FRAC_LEN-1); 
    double interm = diff * factor;
    long long dtol =  interm;
    mant_with_sf = dtol;
    mant_part= mant_with_sf >> (int)sf; 
    regime = (sf_t)sf >> ES;

    //find SREG & REM
    if (regime >= 0) {
        if (regime + 3 < N) {
            SREG = (reg_t)regime + 3;
        } else {
            SREG = (reg_t)N;
        }
    } else {            
        SREG = 2- regime;
    }
    REM=N-SREG;

    if(REM>0){
        if (REM > ES) {
            exponent = sf & ((1 << ES) - 1);
            mantissa= (mantissa_t)mant_part ;
            mantissa.set(FRAC_LEN-1);
        }
        else{
            exponent = sf & ((1 << REM) - 1); 
        }

    }

    result.sign = sign;
    result.regime = regime;
    result.exponent = exponent;
    result.mantissa = mantissa;
    result.isZero = isZero;


    return result;
}
float posit2float(ps_t pos){
	ml_t ml;
	float result=0,mantissa,man_d,u_R,e;

	if(pos.isZero) return result;
	if(pos.regime>=0)	ml=N-(pos.regime+3+ES);
	else	ml=N-2+pos.regime-ES;
	man_d=(float)(1<<ml);
	mantissa=(float)pos.mantissa/man_d;
    
	u_R=hls::pow((float)USEED,(float)pos.regime);
	e=hls::pow(2,pos.exponent);
	result=(float)(u_R*e*(1+mantissa));
	return pos.sign?-result:result;
}
double posit2double(ps_t pos){
	ml_t ml;
    sf_t sf;
	double result=0,mantissa,man_d,u_R,e,m;

	if(pos.isZero) return result;
    if(pos.isInf) return INFINITY;

    mantissa =  (double)pos.mantissa / (double) (1<<(FRAC_LEN-1));
	u_R=hls::pow((double)USEED,pos.regime);
	e=hls::pow(2.0,pos.exponent);
    m= mantissa;
	result=u_R*e*m;

	return pos.sign?-result:result;
}

ps_t int2posit(int in){
	ps_t result;
	regime_t regime=0;
	exponent_t exponent=0;
	mantissa_t mantissa=0;
	bool sign=false, isZero=false;
	reg_t exact,mant_part,SREG,REM;
	sf_t sf;

	double MAX=pow(2,REG_SHIFT*(N-2));//(1<<(REG_SHIFT*(N-2)));

	if(in==0) isZero=true;
	else if(in>=MAX) regime=N-2;
	else if(in<=-MAX) { regime=N-2; sign=1; }
	else{
	//if negative
		if(in<0)  {sign=true; in=-in;}
	// find regime
		sf=(int)log2(in); //note: log2 causes great LUT increase
		exact=1<<sf;
		mant_part=in-exact;
		regime=sf>>ES;
	//find SREG & REM
		SREG=regime>=0?regime+3:2-regime;
		REM=N-SREG;
	// find exponent and mantissa
		if(ES==0){
			if(REM>0) mantissa=mant_part<<(REM-sf);
		}
		else{
			if(REM>=ES) {
				exponent=sf.range(ES-1,0);
				mantissa=mant_part<<(REM-sf-ES);
			}
			else if(REM>0) {
				exponent=sf.range(REM-1,0);
			}
		}
	}
	result.sign=sign;
	result.regime=regime;
	result.exponent=exponent;
	result.mantissa=mantissa;
	result.isZero=isZero;
	return result;
}
ps_t positAdd(ps_t x,ps_t y){

	ps_t result;
	sf_t sf_L,sf_S,sf_x,sf_y,sf_r,diff;
	m_add_t mL,mS,mantissa=0,x_mantissa,y_mantissa;
	bool sign,x_sign,y_sign,ABSxIsGreaterEqual,isZero=0,x_isZero,y_isZero;
    int SA=0;
	regime_t regime,x_regime,y_regime;
	exponent_t exponent=0,x_exponent,y_exponent;
	reg_t SREG,REM;
	//read parameters
	x_sign=x.sign;			 y_sign=y.sign;
	x_regime=x.regime;		 y_regime=y.regime;
	x_exponent=x.exponent;	 y_exponent=y.exponent;
	x_mantissa=x.mantissa;   y_mantissa=y.mantissa;
	x_isZero=x.isZero;		 y_isZero=y.isZero;

	//Set the Flag
	ABSxIsGreaterEqual=	x_regime>y_regime ||(x_regime==y_regime && x_exponent>y_exponent)||
			(x_regime==y_regime && x_exponent==y_exponent &&x_mantissa>=y_mantissa);
	//Sign
	sign=ABSxIsGreaterEqual?x_sign:y_sign;
	//determine sf_x & sf_y &sf_r
	sf_x=((sf_t)x_regime<<ES)+x_exponent;
	sf_y=((sf_t)y_regime<<ES)+y_exponent;
	sf_r=ABSxIsGreaterEqual?sf_x:sf_y;
    sf_L=ABSxIsGreaterEqual?sf_x:sf_y;
    sf_S=ABSxIsGreaterEqual?sf_y:sf_x;
    mL=ABSxIsGreaterEqual?x_mantissa:y_mantissa;
    mS=ABSxIsGreaterEqual?y_mantissa:x_mantissa;


    diff = sf_L -sf_S;
    mS  = mS >>diff;
	if(x_sign ^ y_sign){ //Substraction
        mantissa = mL  - mS;
	}
	else { //Addition  
		mantissa = mL + mS;
	}
    //std::cout<<"mantissa:"<<mantissa<<std::endl;
    SA =  LOD_ADD(mantissa);
    mantissa = mantissa <<SA;
    sf_r = sf_r-SA;

	// set regime & exponent & mantissa
	regime=sf_r>>ES;
	//find SREG & REM
	SREG=regime>=0?regime+3:2-regime;
	REM=N-SREG;
    // Handle exponent and mantissa extraction
    if(REM>0){
        if (REM > ES) {
            exponent = sf_r & ((1 << ES) - 1);
            mantissa= (mantissa_t)mantissa ;
        }
        else{
            exponent = sf_r & ((1 << REM) - 1); 
        }

    }

	if(x_isZero && !y_isZero){ // x=0 ,y!=0 --> r=y
		sign=y_sign;
		regime=y_regime;
		exponent=y_exponent;
		mantissa=y_mantissa;
	}
	else if(!x_isZero && y_isZero){ // y=0 ,x!=0 --> r=x
		sign=x_sign;
		regime=x_regime;
		exponent=x_exponent;
		mantissa=x_mantissa;
	}
    else if((x_isZero && y_isZero)||((x_regime==y_regime && x_exponent==y_exponent && x_mantissa==y_mantissa )&& x.sign!=y.sign)){
		isZero=true;
		regime=0; exponent=0; mantissa=1<<(FRAC_LEN-1);
	}
	result.regime=regime;
	result.exponent=exponent;
	result.mantissa=mantissa;
	result.sign=sign;
	result.isZero=isZero;
	return result;
}

ps_t positDiv2p(ps_t in,int i){
    ps_t result;
    sf_t sf_in, sf_out;
    result = in;
    sf_in=((sf_t)in.regime<<ES)+in.exponent;
    sf_out = sf_in +i;
    if (sf_out.range(REG_LEN+ES,REG_LEN+ES-1)==2){
        result.regime = 2-N;
    }
    else{
        result.regime=sf_out>>ES;
        result.exponent = sf_out & ((1 << ES) - 1);
    }

    return result;

}
ps_t positDiv(ps_t x,ps_t y){

	ps_t result;
	sf_t sf_x,sf_y,sf_r;
	bool sign,x_sign,y_sign,x_isZero,y_isZero,isZero=0;
	regime_t R,regime,x_regime,y_regime;
	exponent_t exponent=0,y_exponent,x_exponent;
	mantissa_t mantissa=0;
	dv_t x_mantissa,y_mantissa,mant=0;
	reg_t SREG,REM;
	int SA=0;
	//read parameters
	x_sign=x.sign;			 y_sign=y.sign;
	x_regime=x.regime;		 y_regime=y.regime;
	x_exponent=x.exponent;	 y_exponent=y.exponent;
	x_mantissa=x.mantissa;   y_mantissa=y.mantissa;
	x_isZero=x.isZero;		 y_isZero=y.isZero;
	//Sign
	sign=x_sign ^ y_sign;
	//isZero
	isZero=x_isZero | y_isZero;
	//determine scaling factors sf_x & sf_y & sf_r
	sf_x=((sf_t)x_regime<<ES)+x_exponent;
	sf_y=((sf_t)y_regime<<ES)+y_exponent;
	sf_r=sf_x-sf_y;
    //Sf Overflow
    if (sf_r.range(REG_LEN+ES,REG_LEN+ES-1)==2){
        regime = 2-N;
    }
    else{

        //multiply mantissa 

        mant = (x_mantissa<<(FRAC_LEN-1))/y_mantissa;

        if(mant[FRAC_LEN-1] == 0){
            mant = (mant<<1);
            //mant[0]= ex_bit;
            sf_r = sf_r-1;
        }
        mantissa = (mantissa_t)mant;   

        //determine regime
        regime=sf_r>>ES;
        //std::cout<<"regime:"<<mant<<std::endl;
    }
    if (sf_r.range(REG_LEN+ES,REG_LEN+ES-1)==2){
        regime = 2-N;
    }
	if(regime>=N-2) regime=N-2;
	else if(regime<2-N) regime=2-N;
	else{
		//find SREG & REM
		SREG=regime>=0?regime+3:2-regime;
		REM=N-SREG;
        if(REM>0){
            if (REM > ES) {
                exponent = sf_r & ((1 << ES) - 1);

            }
            else{
                exponent = sf_r & ((1 << REM) - 1); 
            }

        }
	}
	if(isZero) {regime=0;exponent=0;mantissa=1<<(FRAC_LEN-1);sign=0;}
	//set result
	result.regime=regime;
	result.exponent=exponent;
	result.mantissa=mantissa;
	result.sign=sign;
	result.isZero=isZero;

	return result;
}

ps_t positMul(ps_t x,ps_t y){

	ps_t result;
	sf_t sf_x,sf_y,sf_r,max =((2-N)<<ES);
    
	bool ovf,sign,x_sign,y_sign,x_isZero,y_isZero,isZero=0;
	regime_t R,regime,x_regime,y_regime;
	exponent_t exponent=0,y_exponent,x_exponent;
	mantissa_t mantissa=0;
	mul_t x_mantissa,y_mantissa,mant=0,mul_part=0,mInter=0;
	reg_t SREG,REM;
	int SA=0;
    mantissa.set(FRAC_LEN-1);
	//read parameters
	x_sign=x.sign;			 y_sign=y.sign;
	x_regime=x.regime;		 y_regime=y.regime;
	x_exponent=x.exponent;	 y_exponent=y.exponent;
	x_mantissa=x.mantissa;   y_mantissa=y.mantissa;
	x_isZero=x.isZero;		 y_isZero=y.isZero;
	//Sign
	sign=x_sign ^ y_sign;
	//isZero
	isZero=x_isZero | y_isZero;
	//determine scaling factors sf_x & sf_y & sf_r
	sf_x=((sf_t)x_regime<<ES)+x_exponent;
	sf_y=((sf_t)y_regime<<ES)+y_exponent;
	sf_r=(sf_t)sf_x+sf_y;

    if (sf_r.range(REG_LEN+ES,REG_LEN+ES-1)==2){
        regime = 2-N;
        //std::cout<<"OVF"<<std::endl;
        //std::cout<<"regime:"<<regime<<std::endl;
    }
    else{
        //multiply mantissa 
        mant = (mul_t) x_mantissa * y_mantissa;
        ovf = mant[MUL_LEN-1];
        if (ovf){
            sf_r+=1;
            mantissa = (mantissa_t) mant.range(MUL_LEN-1,MUL_LEN-FRAC_LEN);
        } 
        else
            mantissa = (mantissa_t) mant.range(MUL_LEN-2,MUL_LEN-1-FRAC_LEN);

        //determine regime
        regime=sf_r>>ES;
    }
    //std::cout<<"regime:"<<regime<<std::endl;
	if(regime>=N-2) regime=N-2;
	else if(regime<=2-N) regime=2-N;
	else{
		//find SREG & REM
		SREG=regime>=0?regime+3:2-regime;
		REM=N-SREG;
        if(REM>0){
            if (REM > ES) {
                exponent = sf_r & ((1 << ES) - 1);
                mantissa= (mantissa_t)mantissa ;
            }
            else{
                exponent = sf_r & ((1 << REM) - 1); 
            }

        }
	}
	if(isZero) {regime=0;exponent=0;mantissa=1<<(FRAC_LEN-1);sign=0;}
	//set result
	result.regime=regime;
	result.exponent=exponent;
	result.mantissa=mantissa;
	result.sign=sign;
	result.isZero=isZero;

	return result;
}

ps_t positSub(ps_t x,ps_t y){

	ps_t result;
	y.sign= !(y.sign);
	result=positAdd(x,y);
	return result;
}

///////////////////////////COS////////////////////////////////
double dReduceAngle(double angle, bool &negate) {

    angle = fmod(angle, 2 * PI);
    if (angle > PI) 
        angle -= 2 * PI;
    else if (angle < -PI) 
        angle += 2 * PI;

    negate = false;
    if (angle > PI / 2) {
        angle = PI - angle;
        negate = true;  // Negation für cos(-x) beachten
    } else if (angle < -PI / 2) {
        angle = -PI - angle;
        negate = true;  // Negation für cos(-x) beachten
    }
    return angle;

}

// Taylor-Approximation für cos(x)


double dTailorCos(double x) {
    bool negate;
    double x2,term1,term2,term3,term4;
    x = dReduceAngle(x, negate);  // Reduce the angle to [0, π/2]
    //std::cout<<"dRA: "<<x<<std::endl;
    if(x == PI/2 || x == -PI/2) return 0;
    else if (x== 0) return 1;
    else{
        x2 = x * x;
        term1 = 1.0;
        term2 = x2 / 2.0;

        #if TERMS > 2
            #if APPR_TAILOR ==0
                term3 = x2 * x2 / 24.0;
            #elif APPR_TAILOR ==1
                term3 = x2 * x2 / 32.0;
            #endif
        #endif

        #if TERMS > 3
            #if APPR_TAILOR ==0
                term4 = x2 * term3 / 30.0;
            #elif APPR_TAILOR ==1
                term4 = x2 * term3 / 32.0;
            #endif
        #endif

        #if TERMS == 2
            return negate ? -(term1 - term2) : (term1 - term2);
        #elif TERMS == 3
            return negate ? -(term1 - term2 + term3) : (term1 - term2 + term3);
        #elif TERMS == 4
            return negate ? -(term1 - term2 + term3 - term4) : (term1 - term2 + term3 - term4);
        #else
            return negate ? -(term1 - term2 + term3 - term4) : (term1 - term2 + term3 - term4);
        #endif
    }
 
}

float fReduceAngle(float angle, bool &negate) {
    // Begrenzung auf [-π, π]
    angle = fmod(angle, 2 * PI);
    if (angle > PI) 
        angle -= 2 * PI;
    else if (angle < -PI) 
        angle += 2 * PI;

    // Kosinus-Symmetrie ausnutzen
    negate = false;
    if (angle > PI / 2) {
        angle = PI - angle;
        negate = true;  // Negation für cos(-x) beachten
    } else if (angle < -PI / 2) {
        angle = -PI - angle;
        negate = true;  // Negation für cos(-x) beachten
    }
    
    return angle;
}

float fTailorCos(float x) {
    bool negate;
    float x2,term1,term2,term3,term4;
    x = fReduceAngle(x, negate);  // Reduce the angle to [0, π/2]
    if(x == PI/2 || x == -PI/2) return 0;
    else if (x== 0) return 1;
    else{
        x2 = x * x;
        term1 = 1.0;
        term2 = x2 / 2.0;

        #if TERMS > 2
            #if APPR_TAILOR ==0
                term3 = x2 * x2 / 24.0;
            #elif APPR_TAILOR ==1
                term3 = x2 * x2 / 32.0;
            #endif
        #endif

        #if TERMS > 3
            #if APPR_TAILOR ==0
                term4 = x2 * term3 / 30.0;
            #elif APPR_TAILOR ==1
                term4 = x2 * term3 / 32.0;
            #endif
        #endif

        #if TERMS == 2
            return negate ? -(term1 - term2) : (term1 - term2);
        #elif TERMS == 3
            return negate ? -(term1 - term2 + term3) : (term1 - term2 + term3);
        #elif TERMS == 4
            return negate ? -(term1 - term2 + term3 - term4) : (term1 - term2 + term3 - term4);
        #else
            return negate ? -(term1 - term2 + term3 - term4) : (term1 - term2 + term3 - term4);
        #endif
    }

}
half hReduceAngle(half angle, bool &negate) {
    // Begrenzung auf [-π, π]
    angle = hls::half_fmod(angle, 2 * PI);
    
    if (angle > PI) 
        angle -= 2 * PI;
    else if (angle < -PI) 
        angle += 2 * PI;

    // Kosinus-Symmetrie ausnutzen
    negate = false;
    if (angle > PI / 2) {
        angle = PI - angle;
        negate = true;  // Negation für cos(-x) beachten
    } else if (angle < -PI / 2) {
        angle = -PI - angle;
        negate = true;  // Negation für cos(-x) beachten
    }
    
    return angle;
}
half hTailorCos(half x) {
    bool negate;
    half x2,term1,term2,term3,term4;
    x = hReduceAngle(x, negate);  // Redu ce the angle to [0, π/2]
    if(x == PI/2 || x == -PI/2) return 0;
    else if (x== 0) return 1;
    else{
        x2 = x * x;
        term1 = 1.0;
        term2 = x2 / 2.0;

        #if TERMS > 2
            #if APPR_TAILOR ==0
                term3 = x2 * x2 / 24.0;
            #elif APPR_TAILOR ==1
                term3 = x2 * x2 / 32.0;
            #endif
        #endif

        #if TERMS > 3
            #if APPR_TAILOR ==0
                term4 = x2 * term3 / 30.0;
            #elif APPR_TAILOR ==1
                term4 = x2 * term3 / 32.0;
            #endif
        #endif

        #if TERMS == 2
            return negate ?(half) -(term1 - term2) :(half) (term1 - term2);
        #elif TERMS == 3
            return negate ? (half)-(term1 - term2 + term3) : (half)(term1 - term2 + term3);
        #elif TERMS == 4
            return negate ? (half)-(term1 - term2 + term3 - term4) : (half) (term1 - term2 + term3 - term4);
        #else
            return negate ? (half)-(term1 - term2 + term3 - term4) :(half) (term1 - term2 + term3 - term4);
        #endif
    }

}
#define EPSILON 1e-12  // Adjust this based on precision needs

double safeZero(double value) {
    return (std::abs(value) < EPSILON) ? 0.0 : value;
}
ps_t pReduceAngle(ps_t angle, bool &negate) {
    ps_t m_angle;
    m_angle = positMod(angle,POSIT_2PI);

    if (isGreater_posit(m_angle,POSIT_PI ))
        m_angle = positAdd(m_angle,POSIT_M_2PI);
    else if (isGreater_posit(POSIT_M_PI,m_angle ))  
        m_angle = positAdd(m_angle,POSIT_2PI);

    // Kosinus-Symmetrie ausnutzen
    negate = false;
	if (isGreater_posit(m_angle,POSIT_PI_OVER2 ))	{
		m_angle = positSub(POSIT_PI,m_angle);
        negate = true;  // Negation für cos(-x) beachten
     //   std::cout<<"1 "<<std::endl;
    }else if (isGreater_posit(POSIT_M_PI_OVER2,m_angle )){
		m_angle = positAdd(m_angle,POSIT_PI);
		//m_angle.sign=!m_angle.sign;
        negate = true;  // Negation für cos(-x) beachten
         //   std::cout<<"2 "<<std::endl;
    }
     //   std::cout<<"p negate: "<<negate<<std::endl;
	return m_angle;
}
// Approximate cosine using Taylor series expansion
ps_t positCos(ps_t x) {
    ps_t y, y2, y4, result;
    ps_t term1, term2, term3, t1minust2, term4;
    bool negate;
    term1 = ONE;
    if (x.isZero == 1) return term1;
    y = pReduceAngle(x, negate);
    //std::cout<<"pDG: "<<posit2double(y)<<std::endl;    
    if(y == POSIT_PI_OVER2 || y ==POSIT_M_PI_OVER2) return ZERO;
    else if (y == ZERO) return ONE;
    else{
        y2 = positMul(y, y);
        y4 = positMul(y2, y2);

        term2 = positDiv2p(y2, -1);

        #if TERMS > 2
            #if APPR_TAILOR ==0
                term3 = positDiv(y4, double2posit(24.0));
            #elif APPR_TAILOR ==1
                term3 = positDiv2p(y4, -5);
            #endif
        #endif

        #if TERMS > 3
            #if APPR_TAILOR ==0
                term4 = positDiv(positMul(term3, y2), double2posit(30.0));
            #elif APPR_TAILOR ==1
                term4 = positDiv2p(positMul(term3, y2), -5);
            #endif        
        #endif

        t1minust2 = positSub(term1, term2);

        #if TERMS == 2
            result = t1minust2;
        #elif TERMS == 3
            result = positAdd(t1minust2, term3);
        #elif TERMS == 4
            result = positSub(positAdd(t1minust2, term3), term4);
        #else
            result = positSub(positAdd(t1minust2, term3), term4);
        #endif
        
        return negate ? posit_negate(result) : result;   
    }

}
///////////////////////////SINE////////////////////////////////
double dNAngle(double angle) {
    //std::cout<<"angle: "<<angle<<std::endl;
    angle = fmod(angle, 2 * PI);  // Begrenzung auf [-2π, 2π]
    
    if (angle > PI) 
        angle -= 2 * PI;
    else if (angle < -PI) 
        angle += 2 * PI;

    if (angle > PI / 2) 
        angle = PI - angle;
    else if (angle < -PI / 2) 
        angle = -PI - angle;

    return angle;
}

double dTailorSin(double in) {
    double term1, term2, term3, term4, x;
    x = dNAngle(in);

    term1 = x;
    #if APPR_TAILOR ==0
        term2 = x * x * x / 6.0;
    #elif APPR_TAILOR ==1
        term2 = x * x * x / 8.0;
    #endif
    #if TERMS > 2
        #if APPR_TAILOR ==0
            term3 = term2 * x * x / 20.0;
        #elif APPR_TAILOR ==1
            term3 = term2*x*x/16.0;
        #endif
    #endif
    #if TERMS > 3
        #if APPR_TAILOR ==0
            term4 = term3 * x * x / 42.0;
        #elif APPR_TAILOR ==1
            term4 = term3*x*x/32.0;
        #endif
    #endif

    #if TERMS == 2
        return term1 - term2;
    #elif TERMS == 3
        return term1 - term2 + term3;
    #elif TERMS == 4
        return term1 - term2 + term3 - term4;
    #else
        return term1 - term2 + term3 - term4;
    #endif
}

float fNAngle(float angle) {
    angle = fmod(angle, 2 * PI);  // Begrenzung auf [-2π, 2π]
    
    if (angle > PI) 
        angle -= 2 * PI;
    else if (angle < -PI) 
        angle += 2 * PI;

    if (angle > PI / 2) 
        angle = PI - angle;
    else if (angle < -PI / 2) 
        angle = -PI - angle;

    return angle;
}

float fTailorSin(float in) {
    float term1, term2, term3, term4, x;
    x = fNAngle(in);

    term1 = x;
    #if APPR_TAILOR ==0
        term2 = x * x * x / 6.0;
    #elif APPR_TAILOR ==1
        term2 = x * x * x / 8.0;
    #endif
    #if TERMS > 2
        #if APPR_TAILOR ==0
            term3 = term2 * x * x / 20.0;
        #elif APPR_TAILOR ==1
            term3 = term2*x*x/16.0;
        #endif
    #endif
    #if TERMS > 3
        #if APPR_TAILOR ==0
            term4 = term3 * x * x / 42.0;
        #elif APPR_TAILOR ==1
            term4 = term3*x*x/32.0;
        #endif
    #endif

    #if TERMS == 2
        return term1 - term2;
    #elif TERMS == 3
        return term1 - term2 + term3;
    #elif TERMS == 4
        return term1 - term2 + term3 - term4;
    #else
        return term1 - term2 + term3 - term4;
    #endif
}
half hNAngle(half angle) {
    angle = hls::half_fmod(angle, 2 * PI);  // Begrenzung auf [-2π, 2π]
    
    if (angle > PI) 
        angle -= 2 * PI;
    else if (angle < -PI) 
        angle += 2 * PI;

    if (angle > PI / 2) 
        angle = PI - angle;
    else if (angle < -PI / 2) 
        angle = -PI - angle;

    return angle;
}
half hTailorSin(half in) {
    half term1, term2, term3, term4, x;
    x = hNAngle(in);

    term1 = x;
    #if APPR_TAILOR ==0
        term2 = x * x * x / 6.0;
    #elif APPR_TAILOR ==1
        term2 = x * x * x / 8.0;
    #endif
    #if TERMS > 2
        #if APPR_TAILOR ==0
            term3 = term2 * x * x / 20.0;
        #elif APPR_TAILOR ==1
            term3 = term2*x*x/16.0;
        #endif
    #endif
    #if TERMS > 3
        #if APPR_TAILOR ==0
            term4 = term3 * x * x / 42.0;
        #elif APPR_TAILOR ==1
            term4 = term3*x*x/32.0;
        #endif
    #endif

    #if TERMS == 2
        return term1 - term2;
    #elif TERMS == 3
        return term1 - term2 + term3;
    #elif TERMS == 4
        return term1 - term2 + term3 - term4;
    #else
        return term1 - term2 + term3 - term4;
    #endif
}


ps_t pNAngle(ps_t angle) {
    ps_t m_angle;
    m_angle = positMod(angle,POSIT_2PI);

    if (isGreater_posit(m_angle,POSIT_PI )) 
        m_angle = positAdd(m_angle,POSIT_M_2PI);
	else if (isGreater_posit(POSIT_M_PI,m_angle )) 
		m_angle = positAdd(m_angle,POSIT_2PI);	
	if (isGreater_posit(m_angle,POSIT_PI_OVER2 ))	
		m_angle = positSub(POSIT_PI,m_angle);
	else if (isGreater_posit(POSIT_M_PI_OVER2,m_angle )){
		m_angle = positAdd(m_angle,POSIT_PI);
		m_angle.sign=!m_angle.sign;
        }
	return m_angle;


}


// Approximate sine using Taylor series expansion
ps_t positSin(ps_t x) {
    ps_t y,y2,y3,y5,y7;
    //std::cout<<"x: "<<posit2double(x)<<std::endl;
    y = pNAngle(x);
    //std::cout<<"y: "<<posit2double(y)<<std::endl;
    y2 = positMul(y, y);
    y3 = positMul(y2,y);
    #if TERMS > 2
        y5 = positMul(y2, y3);
    #endif
    #if TERMS > 3
        y7 = positMul(y2, y5);
    #endif
    //2 TERMS
    #if TERMS == 2
        #if APPR_TAILOR ==0
            return positSub(y,positDiv(y3,double2posit(6)));
        #elif APPR_TAILOR ==1
            return positSub(y,positDiv2p(y3,-3));
        #endif

    //3 TERMS
     #elif TERMS == 3
        #if APPR_TAILOR ==0
            return positAdd(positSub(y,positDiv(y3,double2posit(6))),positDiv(y5,double2posit(120)));
        #elif APPR_TAILOR ==1
            return positAdd(positSub(y,positDiv2p(y3,-3)),positDiv2p(y5,-7));
        #endif
    //4 TERMS
    #elif TERMS == 4
        #if APPR_TAILOR ==0
            return positSub(positAdd(positSub(y,positDiv(y3,double2posit(6))),positDiv(y5,double2posit(120))),positDiv(y7, double2posit(5040)));
        #elif APPR_TAILOR ==1
            return positSub(positAdd(positSub(y,positDiv2p(y3,-3)),positDiv2p(y5,-7)),positDiv2p(y7, -12));
        #endif
    #endif
}

void pEuler(ps_t angle, ps_t *result_real, ps_t *result_imag) {
    *result_real = positCos(angle);
    *result_imag = positSin(angle);
}
#if EXACT ==0
void dEuler(double angle, double *result_real, double *result_imag) {
    *result_real = dTailorCos(angle);
    *result_imag = dTailorSin(angle);
}
void fEuler(float angle, float *result_real, float *result_imag) {
    *result_real = fTailorCos(angle);
    *result_imag = fTailorSin(angle);
}
void hEuler(half angle, half *result_real, half *result_imag) {
    *result_real = hTailorCos(angle);
    *result_imag = hTailorSin(angle);
}
#endif
#if EXACT ==1
void dEuler(double angle, double *result_real, double *result_imag) {
    *result_real = hls::cos(angle);
    *result_imag = hls::sin(angle);
}
void fEuler(float angle, float *result_real, float *result_imag) {
    *result_real = hls::cos(angle);
    *result_imag = hls::sin(angle);
}
#endif


void dAccumulateFC(int k, int sampleCount, const double signal[], double& realSum, double& imagSum) {
    double realPart, imagPart, angle = 0.0;
    double deltaTheta = -2.0 * PI * k / sampleCount;

    //addTod_DTH_Array(deltaTheta);
    double signalVal;    

    // Initialize sums

    // Accumulate FFT bin
    for (int n = 0; n < sampleCount; n++) {
        
   /*     // Save current values to files
        realFile << angle<<"   "<<realPart << std::endl;
        imagFile << angle<<"   "<<imagPart << std::endl;*/
        signalVal = signal[n];
        dEuler(angle, &realPart, &imagPart);
        // Accumulate the results
        realSum += signalVal * realPart;
        imagSum += signalVal * imagPart;
        //addTod_Angle_Array(angle);
        angle += deltaTheta;
        
        //addTod_RP_Array(realPart);
		//addTod_IMG_Array(imagPart);
    }
/*    // Close files
    realFile.close();
    imagFile.close();*/
    
}

void fAccumulateFC(int k, int sampleCount, const float signal[], float& realSum, float& imagSum) {
    float realPart, imagPart, angle = 0.0;
    float deltaTheta = -2.0 * PI * k / sampleCount;
    float signalVal;
	//f_addTof_DTH_Array(deltaTheta);
    for (int n = 0; n < sampleCount; n++) {
        signalVal = signal[n];
        fEuler(angle, &realPart, &imagPart);
        realSum += signalVal * realPart;
        imagSum += signalVal * imagPart;
        //f_addTof_Angle_Array(angle);
        angle += deltaTheta;
        //f_addTof_RP_Array(realPart);
		//f_addTof_IMG_Array(imagPart);
    }
}
ps_t calculateKFactor(int k){
	ps_t result;
	log_sf_t mant;
	sf_t sf=0;
	regime_t regime=0;
	exponent_t exponent=0;
	mantissa_t mantissa=0;
	int shiftAmount=0;
	if (k==0) return result;
    result.isZero = false;
	mant = k;

	mantissa.range(FRAC_LEN-2,FRAC_LEN-1-LOG_IN_SIZE_SF) =mant;
	shiftAmount = LOD_FFT(mantissa);
	sf = -shiftAmount;
	regime = (sf_t)sf >> ES;
	exponent = sf & ((1 << ES) - 1);
    //std::cout<<"mantissa: "<<mantissa<<std::endl;
    mantissa = mantissa<<shiftAmount;
	result.regime = regime;
	result.exponent =exponent;
	result.mantissa=mantissa;

	return result;
}

void pAccumulateFC(int k, int sampleCount, const ps_t signal[], ps_t& realSum, ps_t& imagSum) {
    ps_t angle=ZERO, realPart, imagPart;
    ps_t deltaTheta, k_factor;
	k_factor = calculateKFactor(k);
	deltaTheta = positMul(POSIT_M_2PI,k_factor);
    ps_t signalVal;

    for (int n = 0; n < sampleCount; n++) {
        // Access the signal array instead of reading from the stream
        signalVal = signal[n];

        // Compute Euler for posit type

        pEuler(angle, &realPart, &imagPart);
   /*     // Save current values
        realFile << posit2double(angle )<<"   "<<posit2double(realPart) << std::endl;
        imagFile << posit2double(angle )<<"   "<<posit2double(imagPart) << std::endl;*/
        // Accuulate the results
        realSum = positAdd(realSum, positMul(signalVal, realPart));
        imagSum = positAdd(imagSum, positMul(signalVal, imagPart));
        //p_addTop_Angle_Array(posit2double(angle));
        angle = positAdd(angle, deltaTheta);
        
        //p_addTop_RP_Array(posit2double(realPart));
		//p_addTop_IMG_Array(posit2double(imagPart));
    }

}
#if (FFT==0)
// Function for pFFT (using ps_t type)
void pFFT(ps_t signal[], pFFTResult& result) {
    int sampleCount = IN_SIZE;
    ps_t realSum,imagSum;
    for (int k = 0; k < sampleCount; k++) { // For each frequency bin
  //      if (k % 200 == 0)  
   //         std::cout << k << std::endl;

        // Call the accumulation function for each frequency bin
        realSum = ZERO;
        imagSum = ZERO;
        pAccumulateFC(k, sampleCount, signal, realSum, imagSum);

        result.real[k] = realSum;  // Store the real part in the array
        result.imag[k] = imagSum;  // Store the imaginary part in the array
    }
    //deltaThetaPosit.close();
}
#endif
// ====== Bit Reversal ======
void p_bitReverse(ps_t real[], ps_t imag[], int n) {
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (i < j) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
        int m = n / 2;
        while (j >= m && m >= 2) {
            j -= m;
            m /= 2;
        }
        j += m;
    }
}

// ====== In-Place FFT ======
ps_t pMulInt(ps_t p, int i){
    sf_t sf=0;
    ps_t result;
    regime_t regime;
    mantissa_t mantissa;
    exponent_t exponent;
    
    if (i==0) return ZERO;
    else{
        sf= LOD_INT(i); 
        regime = (sf_t)sf >> ES;
        exponent = sf & ((1 << ES) - 1);
        mantissa = i<< (FRAC_LEN-sf-1);
        result.regime = regime;
        result.exponent =exponent;
        result.mantissa=mantissa;  
        result.isZero=0; 
        return positMul(p,result);
    }


}
#if (FFT ==1)
void pFFT(ps_t signal[], pFFTResult& result) {
#pragma HLS INLINE off

    // Step 1: Copy signal into real[] and zero imag[]
    for (int i = 0; i < IN_SIZE; ++i) {
        result.real[i] = signal[i];
        result.imag[i] = ZERO;
    }

    // Step 2: Bit reversal
    p_bitReverse(result.real, result.imag, IN_SIZE);

    // Step 3: FFT main loop
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;  // m = 2^s
        ps_t angleStep = positDiv2p(POSIT_M_2PI,-1*s);
        

        for (int k = 0; k < IN_SIZE; k += m) {
#pragma HLS PIPELINE II=1  // Apply pipeline pragma
            for (int j = 0; j < m / 2; ++j) {
                //ps_t angle = positMul(angleStep , double2posit((double)j));
                ps_t angle = pMulInt(angleStep , j);
                ps_t cosA, sinA;
                pEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                ps_t tReal = positSub(positMul(cosA , result.real[index2]) , positMul(sinA , result.imag[index2]));
                ps_t tImag = positAdd(positMul(sinA , result.real[index2]) , positMul(cosA , result.imag[index2]));

                ps_t uReal = result.real[index1];
                ps_t uImag = result.imag[index1];

                result.real[index1]     = positAdd(uReal , tReal);
                result.imag[index1]     = positAdd(uImag , tImag);
                result.real[index2]     = positSub(uReal , tReal);
                result.imag[index2]     = positSub(uImag , tImag);
            }
        }
    }
}
#endif
// ====== Bit Reversal ======
void d_bitReverse(double real[], double imag[], int n) {
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (i < j) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
        int m = n / 2;
        while (j >= m && m >= 2) {
            j -= m;
            m /= 2;
        }
        j += m;
    }
}
#if (FFT==1)
// ====== In-Place FFT ======
void dFFT(double signal[], dFFTResult& result) {


    // Step 1: Copy signal into real[] and zero imag[]
    for (int i = 0; i < IN_SIZE; ++i) {
        result.real[i] = signal[i];
        result.imag[i] = 0.0f;
    }

    // Step 2: Bit reversal
    d_bitReverse(result.real, result.imag, IN_SIZE);

    // Step 3: FFT main loop
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;  // m = 2^s
        double angleStep = -2.0f * PI / m;

        for (int k = 0; k < IN_SIZE; k += m) {
            for (int j = 0; j < m / 2; ++j) {
                double angle = angleStep * j;
                double cosA, sinA;
                dEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                double tReal = cosA * result.real[index2] - sinA * result.imag[index2];
                double tImag = sinA * result.real[index2] + cosA * result.imag[index2];

                double uReal = result.real[index1];
                double uImag = result.imag[index1];

                result.real[index1]     = uReal + tReal;
                result.imag[index1]     = uImag + tImag;
                result.real[index2]     = uReal - tReal;
                result.imag[index2]     = uImag - tImag;
            }
        }
    }
}
#endif
#if (FFT==0)
// Function for dFFT (using double type)
void dFFT(double signal[], dFFTResult& result) {
    const int sampleCount = IN_SIZE;
    double realSum,imagSum;

    for (int k = 0; k < sampleCount; k++) {
        if (k % 200 == 0)
            std::cout<< k << std::endl;

        realSum = 0.0;
        imagSum = 0.0;

        dAccumulateFC(k, sampleCount, signal, realSum, imagSum);
        //std::cout << "Counter after k=" << k << ": " << globalCounter << std::endl;
        result.real[k] = realSum;  // Store the real part in the array
        result.imag[k] = imagSum;  // Store the imaginary part in the array
    }
    //deltaFile.close();
}
#endif
#if (FFT==0)
// Function for fFFT (using float type)
void fFFT(float signal[], fFFTResult& result) {
    const int sampleCount = IN_SIZE;
    float realSum,imagSum;
    for (int k = 0; k < sampleCount; k++) { // For each frequency bin
        if (k % 200 == 0)  
            std::cout << k << std::endl;

        realSum = 0.0;
        imagSum = 0.0;
        fAccumulateFC(k, sampleCount, signal, realSum, imagSum);

        result.real[k] = realSum;  // Store the real part in the array
        result.imag[k] = imagSum;  // Store the imaginary part in the array
    }
}
#endif
// ====== Bit Reversal ======
void f_bitReverse(float real[], float imag[], int n) {
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (i < j) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
        int m = n / 2;
        while (j >= m && m >= 2) {
            j -= m;
            m /= 2;
        }
        j += m;
    }
}
void h_bitReverse(half real[], half imag[], int n) {
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (i < j) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
        int m = n / 2;
        while (j >= m && m >= 2) {
            j -= m;
            m /= 2;
        }
        j += m;
    }
}
#if (FFT==1)
// ====== In-Place FFT ======
void fFFT(float signal[], fFFTResult& result) {


    // Step 1: Copy signal into real[] and zero imag[]
    for (int i = 0; i < IN_SIZE; ++i) {
        result.real[i] = signal[i];
        result.imag[i] = 0.0f;
    }

    // Step 2: Bit reversal
    f_bitReverse(result.real, result.imag, IN_SIZE);

    // Step 3: FFT main loop
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;  // m = 2^s
        float angleStep = -2.0f * PI / m;

        for (int k = 0; k < IN_SIZE; k += m) {
            for (int j = 0; j < m / 2; ++j) {
                float angle = angleStep * j;
                float cosA, sinA;
                fEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                float tReal = cosA * result.real[index2] - sinA * result.imag[index2];
                float tImag = sinA * result.real[index2] + cosA * result.imag[index2];

                float uReal = result.real[index1];
                float uImag = result.imag[index1];

                result.real[index1]     = uReal + tReal;
                result.imag[index1]     = uImag + tImag;
                result.real[index2]     = uReal - tReal;
                result.imag[index2]     = uImag - tImag;
            }
        }
    }
}
void hFFT(half signal[], hFFTResult& result) {


    // Step 1: Copy signal into real[] and zero imag[]
    for (int i = 0; i < IN_SIZE; ++i) {
        result.real[i] = signal[i];
        result.imag[i] = 0.0f;
    }

    // Step 2: Bit reversal
    h_bitReverse(result.real, result.imag, IN_SIZE);

    // Step 3: FFT main loop
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;  // m = 2^s
        half angleStep = -2.0f * PI / m;

        for (int k = 0; k < IN_SIZE; k += m) {
            for (int j = 0; j < m / 2; ++j) {
                half angle = angleStep * j;
                half cosA, sinA;
                hEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                half tReal = cosA * result.real[index2] - sinA * result.imag[index2];
                half tImag = sinA * result.real[index2] + cosA * result.imag[index2];

                half uReal = result.real[index1];
                half uImag = result.imag[index1];

                result.real[index1]     = uReal + tReal;
                result.imag[index1]     = uImag + tImag;
                result.real[index2]     = uReal - tReal;
                result.imag[index2]     = uImag - tImag;
            }
        }
    }
}
#endif

// Accumulation function for dFFT IFFT (using double arrays)
void dAccumulateFC_IFFT(int k, const double real[], const double imag[], double& realSum, double& imagSum, int sampleCount) {
    realSum = 0.0;
    imagSum = 0.0;
    double realPart, imagPart, angle = 0.0;
    double deltaTheta = 2.0 * PI * k / sampleCount;

    for (int n = 0; n < sampleCount; n++) {
        double signalVal_real = real[n];
        double signalVal_imag = imag[n];

        dEuler(angle, &realPart, &imagPart);

        realSum += signalVal_real * realPart - signalVal_imag * imagPart;
        imagSum += signalVal_real * imagPart + signalVal_imag * realPart;
        angle += deltaTheta;
    }
}
#if (IFFT==0)
// IFFT function for dFFT (using double arrays)
void dIFFT(dFFTResult& result, double signal[], int sampleCount) {
    for (int k = 0; k < sampleCount; k++) {
  /*      if (k % 200 == 0)
            std::cout << k << std::endl;*/

        double realSum = 0.0, imagSum = 0.0;
        dAccumulateFC_IFFT(k, result.real, result.imag, realSum, imagSum, sampleCount);

        signal[k] = realSum / sampleCount;
    }
}
#endif
#if (IFFT==1)
// In-place Inverse FFT
void dIFFT(dFFTResult& result, double signal[], int sampleCount) {



    // Step 2: Bit-reverse
    d_bitReverse(result.real, result.imag, sampleCount);

    // Step 3: Apply FFT with positive twiddle angles
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;
        double angleStep = 2.0f * PI / m;  // <-- positive for IFFT

        for (int k = 0; k < sampleCount; k += m) {
            for (int j = 0; j < m / 2; ++j) {
                double angle = angleStep * j;
                double cosA, sinA;
                dEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                double tReal = cosA * result.real[index2] - sinA * result.imag[index2];
                double tImag = sinA * result.real[index2] + cosA * result.imag[index2];

                double uReal = result.real[index1];
                double uImag = result.imag[index1];

                result.real[index1] = uReal + tReal;
                result.imag[index1] = uImag + tImag;
                result.real[index2] = uReal - tReal;
                result.imag[index2] = uImag - tImag;
            }
        }
    }

    // Step 4: Normalize (divide by IN_SIZE)
    for (int i = 0; i < sampleCount; ++i) {
        signal[i] = result.real[i] / sampleCount;
        // Optionally: discard imag[i] or check if it’s close to 0
    }
}

#endif
// Accumulation function for fFFT IFFT (using float arrays)
void fAccumulateFC_IFFT(int k, const float real[], const float imag[], float& realSum, float& imagSum, int sampleCount) {
    realSum = 0.0f;
    imagSum = 0.0f;
    float realPart, imagPart, angle = 0.0f;
    float deltaTheta = 2.0 * PI * k / sampleCount;

    for (int n = 0; n < sampleCount; n++) {
        float signalVal_real = real[n];
        float signalVal_imag = imag[n];

        fEuler(angle, &realPart, &imagPart);

        realSum += signalVal_real * realPart - signalVal_imag * imagPart;
        imagSum += signalVal_real * imagPart + signalVal_imag * realPart;
        angle += deltaTheta;
    }
}
#if (IFFT==0)
// IFFT function for fFFT (using float arrays)
void fIFFT(fFFTResult& result, float signal[], int sampleCount) {
    for (int k = 0; k < sampleCount; k++) {
 /*       if (k % 200 == 0)
            std::cout << k << std::endl;*/

        float realSum = 0.0f, imagSum = 0.0f;
        fAccumulateFC_IFFT(k, result.real, result.imag, realSum, imagSum, sampleCount);

        signal[k] = realSum / sampleCount;
    }
}
#endif
#if (IFFT==1)
	// In-place Inverse FFT
void fIFFT(fFFTResult& result, float signal[], int sampleCount) {


    // Step 2: Bit-reverse

    f_bitReverse(result.real, result.imag, sampleCount);



    // Step 3: Apply FFT with positive twiddle angles
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;
        float angleStep = 2.0f * PI / m;  // <-- positive for IFFT

        for (int k = 0; k < sampleCount; k += m) {
            for (int j = 0; j < m / 2; ++j) {
                float angle = angleStep * j;
                float cosA, sinA;
                fEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                float tReal = cosA * result.real[index2] - sinA * result.imag[index2];
                float tImag = sinA * result.real[index2] + cosA * result.imag[index2];

                float uReal = result.real[index1];
                float uImag = result.imag[index1];

                result.real[index1] = uReal + tReal;
                result.imag[index1] = uImag + tImag;
                result.real[index2] = uReal - tReal;
                result.imag[index2] = uImag - tImag;
            }
        }
    }

    // Step 4: Normalize (divide by IN_SIZE)
    for (int i = 0; i < sampleCount; ++i) {
        signal[i] = result.real[i] / sampleCount;
        // Optionally: discard imag[i] or check if it’s close to 0
    }
}
void hIFFT(hFFTResult& result, half signal[], int sampleCount) {


    // Step 2: Bit-reverse

    h_bitReverse(result.real, result.imag, sampleCount);



    // Step 3: Apply FFT with positive twiddle angles
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;
        half angleStep = 2.0f * PI / m;  // <-- positive for IFFT

        for (int k = 0; k < sampleCount; k += m) {
            for (int j = 0; j < m / 2; ++j) {
                half angle = angleStep * j;
                half cosA, sinA;
                hEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                half tReal = cosA * result.real[index2] - sinA * result.imag[index2];
                half tImag = sinA * result.real[index2] + cosA * result.imag[index2];

                half uReal = result.real[index1];
                half uImag = result.imag[index1];

                result.real[index1] = uReal + tReal;
                result.imag[index1] = uImag + tImag;
                result.real[index2] = uReal - tReal;
                result.imag[index2] = uImag - tImag;
            }
        }
    }

    // Step 4: Normalize (divide by IN_SIZE)
    for (int i = 0; i < sampleCount; ++i) {
        signal[i] = result.real[i] / sampleCount;
        // Optionally: discard imag[i] or check if it’s close to 0
    }
}
#endif
// Accumulation function for pFFT IFFT (using ps_t arrays)
void pAccumulateFC_IFFT(int k, const ps_t real[], const ps_t imag[], ps_t& realSum, ps_t& imagSum, int sampleCount) {

    ps_t realPart, imagPart, angle = ZERO;
    ps_t deltaTheta, k_factor;
	k_factor = calculateKFactor(k);
	deltaTheta = positMul(POSIT_2PI,k_factor);
    //ps_t deltaTheta = double2posit(2.0 * PI * k / sampleCount);

    for (int n = 0; n < sampleCount; n++) {
        ps_t signalVal_real = real[n];
        ps_t signalVal_imag = imag[n];

        pEuler(angle, &realPart, &imagPart);

        realSum = positAdd(realSum, positSub(positMul(signalVal_real, realPart), positMul(signalVal_imag, imagPart)));
        imagSum = positAdd(imagSum, positAdd(positMul(signalVal_real, imagPart), positMul(signalVal_imag, realPart)));
        angle = positAdd(angle, deltaTheta);
    }

}
#if (IFFT==0)
// IFFT function for pFFT (using ps_t arrays)
void pIFFT(pFFTResult& result, ps_t signal[], int sampleCount) {
    for (int k = 0; k < sampleCount; k++) {
/*        if (k % 200 == 0)
            std::cout << k << std::endl;*/

        ps_t realSum=ZERO;
        ps_t imagSum=ZERO;
        pAccumulateFC_IFFT(k, result.real, result.imag, realSum, imagSum, sampleCount);
        //signal[k] = positDiv(realSum ,double2posit(sampleCount));
		signal[k] = positDiv2p(realSum ,-1*LOG_IN_SIZE_SF);
    }
}
#endif
#if (IFFT==1)
void pIFFT(pFFTResult& result, ps_t signal[], int sampleCount) {


    // Step 2: Bit-reverse
    p_bitReverse(result.real, result.imag, sampleCount);

    // Step 3: Apply FFT with positive twiddle angles
    for (int s = 1; s <= LOG_IN_SIZE_SF; ++s) {
        int m = 1 << s;
        ps_t angleStep = positDiv2p(POSIT_2PI,-1*s);

        for (int k = 0; k < sampleCount; k += m) {
            #pragma HLS PIPELINE II=1  // Apply pipeline pragma
            for (int j = 0; j < m / 2; ++j) {
                ps_t angle = pMulInt(angleStep , j);
                ps_t cosA, sinA;
                pEuler(angle, &cosA, &sinA);

                int index1 = k + j;
                int index2 = k + j + m / 2;

                ps_t tReal = positSub(positMul(cosA , result.real[index2]) , positMul(sinA , result.imag[index2]));
                ps_t tImag = positAdd(positMul(sinA , result.real[index2]) , positMul(cosA , result.imag[index2]));

                ps_t uReal = result.real[index1];
                ps_t uImag = result.imag[index1];

                result.real[index1] = positAdd(uReal , tReal);
                result.imag[index1] = positAdd(uImag , tImag);
                result.real[index2] = positSub(uReal , tReal);
                result.imag[index2] = positSub(uImag , tImag);
            }
        }
    }

    // Step 4: Normalize (divide by IN_SIZE)
    for (int i = 0; i < sampleCount; ++i) {
        signal[i] = positDiv2p(result.real[i] ,-1*LOG_IN_SIZE_SF);
        // Optionally: discard imag[i] or check if it’s close to 0
    }
}	
#endif