/*
  Derived from the Gnu Scientific Library (gsl).
  Most of specfunc/gamma and choose, and dependencies (quite a lot).
*/
/* specfunc/gamma.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* Ilya G. Goldberg: Pulled dependencies for choose and gamma into a single file */
#include "specfunc.h"
#include <math.h>
#include <limits.h>
#include <sys/types.h>

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0)
#define GSL_NEGZERO (-0)


#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02
#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154

#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#define LogRootTwoPi_  0.9189385332046727418
#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif
#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#define GSL_SF_GAMMA_XMAX  171.0
#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

#define GSL_ERROR(reason, gsl_errno)
#define GSL_MAX_DBL(a,b)   fmax(a,b)

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;


#define OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)

#define UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)

#define INTERNAL_OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; return GSL_EOVRFLW; } while(0)

#define INTERNAL_UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; return GSL_EUNDRFLW; } while(0)

#define DOMAIN_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

#define DOMAIN_ERROR_MSG(msg, result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ((msg), GSL_EDOM); } while(0)

#define DOMAIN_ERROR_E10(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; (result)->e10 = 0 ; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

#define OVERFLOW_ERROR_E10(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; (result)->e10 = 0; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
#define UNDERFLOW_ERROR_E10(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; (result)->e10 = 0; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)


#define OVERFLOW_ERROR_2(r1,r2) do { (r1)->val = GSL_POSINF; (r1)->err = GSL_POSINF; (r2)->val = GSL_POSINF ; (r2)->err=GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)

#define UNDERFLOW_ERROR_2(r1,r2) do { (r1)->val = 0.0; (r1)->err = GSL_DBL_MIN; (r2)->val = 0.0 ; (r2)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)

#define DOMAIN_ERROR_2(r1,r2) do { (r1)->val = GSL_NAN; (r1)->err = GSL_NAN;  (r2)->val = GSL_NAN; (r2)->err = GSL_NAN;  GSL_ERROR ("domain error", GSL_EDOM); } while(0)

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))


int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);

/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
  1.00000000000000000000000000000,
  0.083333333333333333333333333333,
 -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
 -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
 -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
 -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
 -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
 -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
 -8.9535174270375468504026113181e-23
};

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

static double apsics_data[16] = {    
  -.0204749044678185,
  -.0101801271534859,
   .0000559718725387,
  -.0000012917176570,
   .0000000572858606,
  -.0000000038213539,
   .0000000003397434,
  -.0000000000374838,
   .0000000000048990,
  -.0000000000007344,
   .0000000000001233,
  -.0000000000000228,
   .0000000000000045,
  -.0000000000000009,
   .0000000000000002,
  -.0000000000000000 
};    
static cheb_series apsi_cs = {
  apsics_data,
  15,
  -1, 1,
  9
};

static double psics_data[23] = {
  -.038057080835217922,
   .491415393029387130, 
  -.056815747821244730,
   .008357821225914313,
  -.001333232857994342,
   .000220313287069308,
  -.000037040238178456,
   .000006283793654854,
  -.000001071263908506,
   .000000183128394654,
  -.000000031353509361,
   .000000005372808776,
  -.000000000921168141,
   .000000000157981265,
  -.000000000027098646,
   .000000000004648722,
  -.000000000000797527,
   .000000000000136827,
  -.000000000000023475,
   .000000000000004027,
  -.000000000000000691,
   .000000000000000118,
  -.000000000000000020
};
static cheb_series psi_cs = {
  psics_data,
  22,
  -1, 1,
  17
};


#define PSI_TABLE_NMAX 100
static double psi_table[PSI_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(0) */
 -M_EULER,                          /* psi(1) */
  0.42278433509846713939348790992,  /* ...    */
  0.92278433509846713939348790992,
  1.25611766843180047272682124325,
  1.50611766843180047272682124325,
  1.70611766843180047272682124325,
  1.87278433509846713939348790992,
  2.01564147795560999653634505277,
  2.14064147795560999653634505277,
  2.25175258906672110764745616389,
  2.35175258906672110764745616389,
  2.44266167997581201673836525479,
  2.52599501330914535007169858813,
  2.60291809023222227314862166505,
  2.67434666166079370172005023648,
  2.74101332832746036838671690315,
  2.80351332832746036838671690315,
  2.86233685773922507426906984432,
  2.91789241329478062982462539988,
  2.97052399224214905087725697883,
  3.02052399224214905087725697883,
  3.06814303986119666992487602645,
  3.11359758531574212447033057190,
  3.15707584618530734186163491973,
  3.1987425128519740085283015864,
  3.2387425128519740085283015864,
  3.2772040513135124700667631249,
  3.3142410883505495071038001619,
  3.3499553740648352213895144476,
  3.3844381326855248765619282407,
  3.4177714660188582098952615740,
  3.4500295305349872421533260902,
  3.4812795305349872421533260902,
  3.5115825608380175451836291205,
  3.5409943255438998981248055911,
  3.5695657541153284695533770196,
  3.5973435318931062473311547974,
  3.6243705589201332743581818244,
  3.6506863483938174848844976139,
  3.6763273740348431259101386396,
  3.7013273740348431259101386396,
  3.7257176179372821503003825420,
  3.7495271417468059598241920658,
  3.7727829557002943319172153216,
  3.7955102284275670591899425943,
  3.8177324506497892814121648166,
  3.8394715810845718901078169905,
  3.8607481768292527411716467777,
  3.8815815101625860745049801110,
  3.9019896734278921969539597029,
  3.9219896734278921969539597029,
  3.9415975165651470989147440166,
  3.9608282857959163296839747858,
  3.9796962103242182164764276160,
  3.9982147288427367349949461345,
  4.0163965470245549168131279527,
  4.0342536898816977739559850956,
  4.0517975495308205809735289552,
  4.0690389288411654085597358518,
  4.0859880813835382899156680552,
  4.1026547480502049565823347218,
  4.1190481906731557762544658694,
  4.1351772229312202923834981274,
  4.1510502388042361653993711433,
  4.1666752388042361653993711433,
  4.1820598541888515500147557587,
  4.1972113693403667015299072739,
  4.2121367424746950597388624977,
  4.2268426248276362362094507330,
  4.2413353784508246420065521823,
  4.2556210927365389277208378966,
  4.2697055997787924488475984600,
  4.2835944886676813377364873489,
  4.2972931188046676391063503626,
  4.3108066323181811526198638761,
  4.3241399656515144859531972094,
  4.3372978603883565912163551041,
  4.3502848733753695782293421171,
  4.3631053861958823987421626300,
  4.3757636140439836645649474401,
  4.3882636140439836645649474401,
  4.4006092930563293435772931191,
  4.4128044150075488557724150703,
  4.4248526077786331931218126607,
  4.4367573696833950978837174226,
  4.4485220755657480390601880108,
  4.4601499825424922251066996387,
  4.4716442354160554434975042364,
  4.4830078717796918071338678728,
  4.4942438268358715824147667492,
  4.5053549379469826935258778603,
  4.5163439489359936825368668713,
  4.5272135141533849868846929582,
  4.5379662023254279976373811303,
  4.5486045001977684231692960239,
  4.5591308159872421073798223397,
  4.5695474826539087740464890064,
  4.5798567610044242379640147796,
  4.5900608426370772991885045755,
  4.6001618527380874001986055856
};

#define PSI_1_TABLE_NMAX 100
static double psi_1_table[PSI_1_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(1,0) */
  M_PI*M_PI/6.0,                    /* psi(1,1) */
  0.644934066848226436472415,       /* ...      */
  0.394934066848226436472415,
  0.2838229557371153253613041,
  0.2213229557371153253613041,
  0.1813229557371153253613041,
  0.1535451779593375475835263,
  0.1331370146940314251345467,
  0.1175120146940314251345467,
  0.1051663356816857461222010,
  0.0951663356816857461222010,
  0.0869018728717683907503002,
  0.0799574284273239463058557,
  0.0740402686640103368384001,
  0.0689382278476838062261552,
  0.0644937834032393617817108,
  0.0605875334032393617817108,
  0.0571273257907826143768665,
  0.0540409060376961946237801,
  0.0512708229352031198315363,
  0.0487708229352031198315363,
  0.0465032492390579951149830,
  0.0444371335365786562720078,
  0.0425467743683366902984728,
  0.0408106632572255791873617,
  0.0392106632572255791873617,
  0.0377313733163971768204978,
  0.0363596312039143235969038,
  0.0350841209998326909438426,
  0.0338950603577399442137594,
  0.0327839492466288331026483,
  0.0317433665203020901265817,
  0.03076680402030209012658168,
  0.02984853037475571730748159,
  0.02898347847164153045627052,
  0.02816715194102928555831133,
  0.02739554700275768062003973,
  0.02666508681283803124093089,
  0.02597256603721476254286995,
  0.02531510384129102815759710,
  0.02469010384129102815759710,
  0.02409521984367056414807896,
  0.02352832641963428296894063,
  0.02298749353699501850166102,
  0.02247096461137518379091722,
  0.02197713745088135663042339,
  0.02150454765882086513703965,
  0.02105185413233829383780923,
  0.02061782635456051606003145,
  0.02020133322669712580597065,
  0.01980133322669712580597065,
  0.01941686571420193164987683,
  0.01904704322899483105816086,
  0.01869104465298913508094477,
  0.01834810912486842177504628,
  0.01801753061247172756017024,
  0.01769865306145131939690494,
  0.01739086605006319997554452,
  0.01709360088954001329302371,
  0.01680632711763538818529605,
  0.01652854933985761040751827,
  0.01625980437882562975715546,
  0.01599965869724394401313881,
  0.01574770606433893015574400,
  0.01550356543933893015574400,
  0.01526687904880638577704578,
  0.01503731063741979257227076,
  0.01481454387422086185273411,
  0.01459828089844231513993134,
  0.01438824099085987447620523,
  0.01418415935820681325171544,
  0.01398578601958352422176106,
  0.01379288478501562298719316,
  0.01360523231738567365335942,
  0.01342261726990576130858221,
  0.01324483949212798353080444,
  0.01307170929822216635628920,
  0.01290304679189732236910755,
  0.01273868124291638877278934,
  0.01257845051066194236996928,
  0.01242220051066194236996928,
  0.01226978472038606978956995,
  0.01212106372098095378719041,
  0.01197590477193174490346273,
  0.01183418141592267460867815,
  0.01169577311142440471248438,
  0.01156056489076458859566448,
  0.01142844704164317229232189,
  0.01129931481023821361463594,
  0.01117306812421372175754719,
  0.01104961133409026496742374,
  0.01092885297157366069257770,
  0.01081070552355853781923177,
  0.01069508522063334415522437,
  0.01058191183901270133041676,
  0.01047110851491297833872701,
  0.01036260157046853389428257,
  0.01025632035036012704977199,  /* ...        */
  0.01015219706839427948625679,  /* psi(1,99)  */
  0.01005016666333357139524567   /* psi(1,100) */
};


/* coefficients for gamma=7, kmax=8  Lanczos method */
static double lanczos_7_c[9] = {
  0.99999999999980993227684700473478,
  676.520368121885098567009190444019,
 -1259.13921672240287047156078755283,
  771.3234287776530788486528258894,
 -176.61502916214059906584551354,
  12.507343278686904814458936853,
 -0.13857109526572011689554707,
  9.984369578019570859563e-6,
  1.50563273514931155834e-7
};
/* Chebyshev expansion for log(gamma(x)/gamma(8))
 * 5 < x < 10
 * -1 < t < 1
 */
static double gamma_5_10_data[24] = {
 -1.5285594096661578881275075214,
  4.8259152300595906319768555035,
  0.2277712320977614992970601978,
 -0.0138867665685617873604917300,
  0.0012704876495201082588139723,
 -0.0001393841240254993658962470,
  0.0000169709242992322702260663,
 -2.2108528820210580075775889168e-06,
  3.0196602854202309805163918716e-07,
 -4.2705675000079118380587357358e-08,
  6.2026423818051402794663551945e-09,
 -9.1993973208880910416311405656e-10,
  1.3875551258028145778301211638e-10,
 -2.1218861491906788718519522978e-11,
  3.2821736040381439555133562600e-12,
 -5.1260001009953791220611135264e-13,
  8.0713532554874636696982146610e-14,
 -1.2798522376569209083811628061e-14,
  2.0417711600852502310258808643e-15,
 -3.2745239502992355776882614137e-16,
  5.2759418422036579482120897453e-17,
 -8.5354147151695233960425725513e-18,
  1.3858639703888078291599886143e-18,
 -2.2574398807738626571560124396e-19
};
static const cheb_series gamma_5_10_cs = {
  gamma_5_10_data,
  23,
  -1, 1,
  11
};


#define GSL_SF_FACT_NMAX 170
static struct {int n; double f; long i; } fact_table[GSL_SF_FACT_NMAX + 1] = {
    { 0,  1.0,     1L     },
    { 1,  1.0,     1L     },
    { 2,  2.0,     2L     },
    { 3,  6.0,     6L     },
    { 4,  24.0,    24L    },
    { 5,  120.0,   120L   },
    { 6,  720.0,   720L   },
    { 7,  5040.0,  5040L  },
    { 8,  40320.0, 40320L },

    { 9,  362880.0,     362880L    },
    { 10, 3628800.0,    3628800L   },
    { 11, 39916800.0,   39916800L  },
    { 12, 479001600.0,  479001600L },

    { 13, 6227020800.0,                               0 },
    { 14, 87178291200.0,                              0 },
    { 15, 1307674368000.0,                            0 },
    { 16, 20922789888000.0,                           0 },
    { 17, 355687428096000.0,                          0 },
    { 18, 6402373705728000.0,                         0 },
    { 19, 121645100408832000.0,                       0 },
    { 20, 2432902008176640000.0,                      0 },
    { 21, 51090942171709440000.0,                     0 },
    { 22, 1124000727777607680000.0,                   0 },
    { 23, 25852016738884976640000.0,                  0 },
    { 24, 620448401733239439360000.0,                 0 },
    { 25, 15511210043330985984000000.0,               0 },
    { 26, 403291461126605635584000000.0,              0 },
    { 27, 10888869450418352160768000000.0,            0 },
    { 28, 304888344611713860501504000000.0,           0 },
    { 29, 8841761993739701954543616000000.0,          0 },
    { 30, 265252859812191058636308480000000.0,        0 },
    { 31, 8222838654177922817725562880000000.0,       0 },
    { 32, 263130836933693530167218012160000000.0,     0 },
    { 33, 8683317618811886495518194401280000000.0,    0 },
    { 34, 2.95232799039604140847618609644e38,  0 },
    { 35, 1.03331479663861449296666513375e40,  0 },
    { 36, 3.71993326789901217467999448151e41,  0 },
    { 37, 1.37637530912263450463159795816e43,  0 },
    { 38, 5.23022617466601111760007224100e44,  0 },
    { 39, 2.03978820811974433586402817399e46,  0 },
    { 40, 8.15915283247897734345611269600e47,  0 },
    { 41, 3.34525266131638071081700620534e49,  0 },
    { 42, 1.40500611775287989854314260624e51,  0 },
    { 43, 6.04152630633738356373551320685e52,  0 },
    { 44, 2.65827157478844876804362581101e54,  0 },
    { 45, 1.19622220865480194561963161496e56,  0 },
    { 46, 5.50262215981208894985030542880e57,  0 },
    { 47, 2.58623241511168180642964355154e59,  0 },
    { 48, 1.24139155925360726708622890474e61,  0 },
    { 49, 6.08281864034267560872252163321e62,  0 },
    { 50, 3.04140932017133780436126081661e64,  0 },
    { 51, 1.55111875328738228022424301647e66,  0 },
    { 52, 8.06581751709438785716606368564e67,  0 },
    { 53, 4.27488328406002556429801375339e69,  0 },
    { 54, 2.30843697339241380472092742683e71,  0 },
    { 55, 1.26964033536582759259651008476e73,  0 },
    { 56, 7.10998587804863451854045647464e74,  0 },
    { 57, 4.05269195048772167556806019054e76,  0 },
    { 58, 2.35056133128287857182947491052e78,  0 },
    { 59, 1.38683118545689835737939019720e80,  0 },
    { 60, 8.32098711274139014427634118320e81,  0 },
    { 61, 5.07580213877224798800856812177e83,  0 },
    { 62, 3.14699732603879375256531223550e85,  0 },
    { 63, 1.982608315404440064116146708360e87,  0 },
    { 64, 1.268869321858841641034333893350e89,  0 },
    { 65, 8.247650592082470666723170306800e90,  0 },
    { 66, 5.443449390774430640037292402480e92,  0 },
    { 67, 3.647111091818868528824985909660e94,  0 },
    { 68, 2.480035542436830599600990418570e96,  0 },
    { 69, 1.711224524281413113724683388810e98,  0 },
    { 70, 1.197857166996989179607278372170e100,  0 },
    { 71, 8.504785885678623175211676442400e101,  0 },
    { 72, 6.123445837688608686152407038530e103,  0 },
    { 73, 4.470115461512684340891257138130e105,  0 },
    { 74, 3.307885441519386412259530282210e107,  0 },
    { 75, 2.480914081139539809194647711660e109,  0 },
    { 76, 1.885494701666050254987932260860e111,  0 },
    { 77, 1.451830920282858696340707840860e113,  0 },
    { 78, 1.132428117820629783145752115870e115,  0 },
    { 79, 8.946182130782975286851441715400e116,  0 },
    { 80, 7.156945704626380229481153372320e118,  0 },
    { 81, 5.797126020747367985879734231580e120,  0 },
    { 82, 4.753643337012841748421382069890e122,  0 },
    { 83, 3.945523969720658651189747118010e124,  0 },
    { 84, 3.314240134565353266999387579130e126,  0 },
    { 85, 2.817104114380550276949479442260e128,  0 },
    { 86, 2.422709538367273238176552320340e130,  0 },
    { 87, 2.107757298379527717213600518700e132,  0 },
    { 88, 1.854826422573984391147968456460e134,  0 },
    { 89, 1.650795516090846108121691926250e136,  0 },
    { 90, 1.485715964481761497309522733620e138,  0 },
    { 91, 1.352001527678402962551665687590e140,  0 },
    { 92, 1.243841405464130725547532432590e142,  0 },
    { 93, 1.156772507081641574759205162310e144,  0 },
    { 94, 1.087366156656743080273652852570e146,  0 },
    { 95, 1.032997848823905926259970209940e148,  0 },
    { 96, 9.916779348709496892095714015400e149,  0 },
    { 97, 9.619275968248211985332842594960e151,  0 },
    { 98, 9.426890448883247745626185743100e153,  0 },
    { 99, 9.332621544394415268169923885600e155,  0 },
    { 100, 9.33262154439441526816992388563e157,  0 },
    { 101, 9.42594775983835942085162312450e159,  0 },
    { 102, 9.61446671503512660926865558700e161,  0 },
    { 103, 9.90290071648618040754671525458e163,  0 },
    { 104, 1.02990167451456276238485838648e166,  0 },
    { 105, 1.08139675824029090050410130580e168,  0 },
    { 106, 1.146280563734708354534347384148e170,  0 },
    { 107, 1.226520203196137939351751701040e172,  0 },
    { 108, 1.324641819451828974499891837120e174,  0 },
    { 109, 1.443859583202493582204882102460e176,  0 },
    { 110, 1.588245541522742940425370312710e178,  0 },
    { 111, 1.762952551090244663872161047110e180,  0 },
    { 112, 1.974506857221074023536820372760e182,  0 },
    { 113, 2.231192748659813646596607021220e184,  0 },
    { 114, 2.543559733472187557120132004190e186,  0 },
    { 115, 2.925093693493015690688151804820e188,  0 },
    { 116, 3.393108684451898201198256093590e190,  0 },
    { 117, 3.96993716080872089540195962950e192,  0 },
    { 118, 4.68452584975429065657431236281e194,  0 },
    { 119, 5.57458576120760588132343171174e196,  0 },
    { 120, 6.68950291344912705758811805409e198,  0 },
    { 121, 8.09429852527344373968162284545e200,  0 },
    { 122, 9.87504420083360136241157987140e202,  0 },
    { 123, 1.21463043670253296757662432419e205,  0 },
    { 124, 1.50614174151114087979501416199e207,  0 },
    { 125, 1.88267717688892609974376770249e209,  0 },
    { 126, 2.37217324288004688567714730514e211,  0 },
    { 127, 3.01266001845765954480997707753e213,  0 },
    { 128, 3.85620482362580421735677065923e215,  0 },
    { 129, 4.97450422247728744039023415041e217,  0 },
    { 130, 6.46685548922047367250730439554e219,  0 },
    { 131, 8.47158069087882051098456875820e221,  0 },
    { 132, 1.11824865119600430744996307608e224,  0 },
    { 133, 1.48727070609068572890845089118e226,  0 },
    { 134, 1.99294274616151887673732419418e228,  0 },
    { 135, 2.69047270731805048359538766215e230,  0 },
    { 136, 3.65904288195254865768972722052e232,  0 },
    { 137, 5.01288874827499166103492629211e234,  0 },
    { 138, 6.91778647261948849222819828311e236,  0 },
    { 139, 9.61572319694108900419719561353e238,  0 },
    { 140, 1.34620124757175246058760738589e241,  0 },
    { 141, 1.89814375907617096942852641411e243,  0 },
    { 142, 2.69536413788816277658850750804e245,  0 },
    { 143, 3.85437071718007277052156573649e247,  0 },
    { 144, 5.55029383273930478955105466055e249,  0 },
    { 145, 8.04792605747199194484902925780e251,  0 },
    { 146, 1.17499720439091082394795827164e254,  0 },
    { 147, 1.72724589045463891120349865931e256,  0 },
    { 148, 2.55632391787286558858117801578e258,  0 },
    { 149, 3.80892263763056972698595524351e260,  0 },
    { 150, 5.71338395644585459047893286526e262,  0 },
    { 151, 8.62720977423324043162318862650e264,  0 },
    { 152, 1.31133588568345254560672467123e267,  0 },
    { 153, 2.00634390509568239477828874699e269,  0 },
    { 154, 3.08976961384735088795856467036e271,  0 },
    { 155, 4.78914290146339387633577523906e273,  0 },
    { 156, 7.47106292628289444708380937294e275,  0 },
    { 157, 1.17295687942641442819215807155e278,  0 },
    { 158, 1.85327186949373479654360975305e280,  0 },
    { 159, 2.94670227249503832650433950735e282,  0 },
    { 160, 4.71472363599206132240694321176e284,  0 },
    { 161, 7.59070505394721872907517857094e286,  0 },
    { 162, 1.22969421873944943411017892849e289,  0 },
    { 163, 2.00440157654530257759959165344e291,  0 },
    { 164, 3.28721858553429622726333031164e293,  0 },
    { 165, 5.42391066613158877498449501421e295,  0 },
    { 166, 9.00369170577843736647426172359e297,  0 },
    { 167, 1.50361651486499904020120170784e300,  0 },
    { 168, 2.52607574497319838753801886917e302,  0 },
    { 169, 4.26906800900470527493925188890e304,  0 },
    { 170, 7.25741561530799896739672821113e306,  0 },

    /*
    { 171, 1.24101807021766782342484052410e309,  0 },
    { 172, 2.13455108077438865629072570146e311,  0 },
    { 173, 3.69277336973969237538295546352e313,  0 },
    { 174, 6.42542566334706473316634250653e315,  0 },
    { 175, 1.12444949108573632830410993864e318,  0 },
    { 176, 1.97903110431089593781523349201e320,  0 },
    { 177, 3.50288505463028580993296328086e322,  0 },
    { 178, 6.23513539724190874168067463993e324,  0 },
    { 179, 1.11608923610630166476084076055e327,  0 },
    { 180, 2.00896062499134299656951336898e329,  0 },
    { 181, 3.63621873123433082379081919786e331,  0 },
    { 182, 6.61791809084648209929929094011e333,  0 },
    { 183, 1.21107901062490622417177024204e336,  0 },
    { 184, 2.22838537954982745247605724535e338,  0 },
    { 185, 4.12251295216718078708070590390e340,  0 },
    { 186, 7.66787409103095626397011298130e342,  0 },
    { 187, 1.43389245502278882136241112750e345,  0 },
    { 188, 2.69571781544284298416133291969e347,  0 },
    { 189, 5.09490667118697324006491921822e349,  0 },
    { 190, 9.68032267525524915612334651460e351,  0 },
    { 191, 1.84894163097375258881955918429e354,  0 },
    { 192, 3.54996793146960497053355363384e356,  0 },
    { 193, 6.85143810773633759312975851330e358,  0 },
    { 194, 1.32917899290084949306717315158e361,  0 },
    { 195, 2.59189903615665651148098764559e363,  0 },
    { 196, 5.08012211086704676250273578535e365,  0 },
    { 197, 1.00078405584080821221303894971e368,  0 },
    { 198, 1.98155243056480026018181712043e370,  0 },
    { 199, 3.94328933682395251776181606966e372,  0 },
    { 200, 7.88657867364790503552363213932e374,  0 }
    */
};





static inline int
cheb_eval_e(const cheb_series * cs,
            const double x,
            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  { 
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}

int gsl_sf_exp_mult_err_e(const double x, const double dx,
                             const double y, const double dy,
                             gsl_sf_result * result)
{
  const double ay  = fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = fabs(dy * exp(x));
    return GSL_SUCCESS;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
          && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    double ex = exp(x);
    result->val  = y * ex;
    result->err  = ex * (fabs(dy) + fabs(y*dx));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double ly  = log(ay);
    const double lnr = x + ly;

    if(lnr > GSL_LOG_DBL_MAX - 0.01) {
      OVERFLOW_ERROR(result);
    }
    else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
      UNDERFLOW_ERROR(result);
    }
    else {
      const double sy  = GSL_SIGN(y);
      const double M   = floor(x);
      const double N   = floor(ly);
      const double a   = x  - M;
      const double b   = ly - N;
      const double eMN = exp(M+N);
      const double eab = exp(a+b);
      result->val  = sy * eMN * eab;
      result->err  = eMN * eab * 2.0*GSL_DBL_EPSILON;
      result->err += eMN * eab * fabs(dy/y);
      result->err += eMN * eab * fabs(dx);
      return GSL_SUCCESS;
    }
  }
  return GSL_SUCCESS;
}


int
gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result)
{
  const double adx = fabs(dx);

  /* CHECK_POINTER(result) */

  if(x + adx > GSL_LOG_DBL_MAX) {
    OVERFLOW_ERROR(result);
  }
  else if(x - adx < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else {
    const double ex  = exp(x);
    const double edx = exp(adx);
    result->val  = ex;
    result->err  = ex * GSL_MAX_DBL(GSL_DBL_EPSILON, edx - 1.0/edx);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  return GSL_SUCCESS;
}


int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s <= 1.0 || q <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    const double max_bits = 54.0;
    const double ln_term0 = -s * log(q);  

    if(ln_term0 < GSL_LOG_DBL_MIN + 1.0) {
      UNDERFLOW_ERROR(result);
    }
    else if(ln_term0 > GSL_LOG_DBL_MAX - 1.0) {
      OVERFLOW_ERROR (result);
    }
    else if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25)) {
      result->val = pow(q, -s);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(s > 0.5*max_bits && q < 1.0) {
      const double p1 = pow(q, -s);
      const double p2 = pow(q/(1.0+q), s);
      const double p3 = pow(q/(2.0+q), s);
      result->val = p1 * (1.0 + p2 + p3);
      result->err = GSL_DBL_EPSILON * (0.5*s + 2.0) * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Euler-Maclaurin summation formula 
       * [Moshier, p. 400, with several typo corrections]
       */
      const int jmax = 12;
      const int kmax = 10;
      int j, k;
      const double pmax  = pow(kmax + q, -s);
      double scp = s;
      double pcp = pmax / (kmax + q);
      double ans = pmax*((kmax+q)/(s-1.0) + 0.5);

      for(k=0; k<kmax; k++) {
        ans += pow(k + q, -s);
      }

      for(j=0; j<=jmax; j++) {
        double delta = hzeta_c[j+1] * scp * pcp;
        ans += delta;
        if(fabs(delta/ans) < 0.5*GSL_DBL_EPSILON) break;
        scp *= (s+2*j+1)*(s+2*j+2);
        pcp /= (kmax + q)*(kmax + q);
      }

      result->val = ans;
      result->err = 2.0 * (jmax + 1.0) * GSL_DBL_EPSILON * fabs(ans);
      return GSL_SUCCESS;
    }
  }
  return GSL_SUCCESS;
}


/* digamma for x both positive and negative; we do both
 * cases here because of the way we use even/odd parts
 * of the function
 */
static int
psi_x(const double x, gsl_sf_result * result)
{
  const double y = fabs(x);

  if(x == 0.0 || x == -1.0 || x == -2.0) {
    DOMAIN_ERROR(result);
  }
  else if(y >= 2.0) {
    const double t = 8.0/(y*y)-1.0;
    gsl_sf_result result_c;
    cheb_eval_e(&apsi_cs, t, &result_c);
    if(x < 0.0) {
      const double s = sin(M_PI*x);
      const double c = cos(M_PI*x);
      if(fabs(s) < 2.0*GSL_SQRT_DBL_MIN) {
        DOMAIN_ERROR(result);
      }
      else {
        result->val  = log(y) - 0.5/x + result_c.val - M_PI * c/s;
        result->err  = M_PI*fabs(x)*GSL_DBL_EPSILON/(s*s);
        result->err += result_c.err;
        result->err += GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
    }
    else {
      result->val  = log(y) - 0.5/x + result_c.val;
      result->err  = result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }
  else { /* -2 < x < 2 */
    gsl_sf_result result_c;

    if(x < -1.0) { /* x = -2 + v */
      const double v  = x + 2.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/(x+1.0);
      const double t3 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);
      
      result->val  = -(t1 + t2 + t3) + result_c.val;
      result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)) + fabs(x/(t3*t3)));
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < 0.0) { /* x = -1 + v */
      const double v  = x + 1.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);
      
      result->val  = -(t1 + t2) + result_c.val;
      result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)));
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < 1.0) { /* x = v */
      const double t1 = 1.0/x;
      cheb_eval_e(&psi_cs, 2.0*x-1.0, &result_c);
      
      result->val  = -t1 + result_c.val;
      result->err  = GSL_DBL_EPSILON * t1;
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else { /* x = 1 + v */
      const double v = x - 1.0;
      return cheb_eval_e(&psi_cs, 2.0*v-1.0, result);
    }
  }
  return GSL_SUCCESS;
}


int gsl_sf_psi_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  return psi_x(x, result);
}


/* generic polygamma; assumes n >= 0 and x > 0
 */
static int
psi_n_xg0(const int n, const double x, gsl_sf_result * result)
{
  if(n == 0) {
    return gsl_sf_psi_e(x, result);
  }
  else {
    /* Abramowitz + Stegun 6.4.10 */
    gsl_sf_result ln_nf;
    gsl_sf_result hzeta;
    int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
    int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
    int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
                                           hzeta.val, hzeta.err,
                                           result);
    if(GSL_IS_EVEN(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
  }
}


int gsl_sf_psi_1_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0 || x == -1.0 || x == -2.0) {
    DOMAIN_ERROR(result);
  }
  else if(x > 0.0)
  {
    return psi_n_xg0(1, x, result);
  }
  else if(x > -5.0)
  {
    /* Abramowitz + Stegun 6.4.6 */
    int M = (int)-floor(x);
    double fx = x + M;
    double sum = 0.0;
    int m;

    if(fx == 0.0)
      DOMAIN_ERROR(result);

    for(m = 0; m < M; ++m)
      sum += 1.0/((x+m)*(x+m));

    {
      int stat_psi = psi_n_xg0(1, fx, result);
      result->val += sum;
      result->err += M * GSL_DBL_EPSILON * sum;
      return stat_psi;
    }
  }
  else
  {
    /* Abramowitz + Stegun 6.4.7 */
    const double sin_px = sin(M_PI * x);
    const double d = M_PI*M_PI/(sin_px*sin_px);
    gsl_sf_result r;
    int stat_psi = psi_n_xg0(1, 1.0-x, &r);
    result->val = d - r.val;
    result->err = r.err + 2.0*GSL_DBL_EPSILON*d;
    return stat_psi;
  }
  return GSL_SUCCESS;
}


int gsl_sf_psi_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= 0) {
    DOMAIN_ERROR(result);
  }
  else if(n <= PSI_TABLE_NMAX) {
    result->val = psi_table[n];
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    /* Abramowitz+Stegun 6.3.18 */
    const double c2 = -1.0/12.0;
    const double c3 =  1.0/120.0;
    const double c4 = -1.0/252.0;
    const double c5 =  1.0/240.0;
    const double ni2 = (1.0/n)*(1.0/n);
    const double ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2*c5)));
    result->val  = log(n) - 0.5/n + ser;
    result->err  = GSL_DBL_EPSILON * (fabs(log(n)) + fabs(0.5/n) + fabs(ser));
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  return GSL_SUCCESS;
}

int gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n == 0)
  {
    return gsl_sf_psi_e(x, result);
  }
  else if(n == 1)
  {
    return gsl_sf_psi_1_e(x, result);
  }
  else if(n < 0 || x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    gsl_sf_result ln_nf;
    gsl_sf_result hzeta;
    int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
    int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
    int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
                                           hzeta.val, hzeta.err,
                                           result);
    if(GSL_IS_EVEN(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
  }
  return GSL_SUCCESS;
}


int gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  if(n <= 0) {
    DOMAIN_ERROR(result);
  }
  else if(n <= PSI_1_TABLE_NMAX) {
    result->val = psi_1_table[n];
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    /* Abramowitz+Stegun 6.4.12
     * double-precision for n > 100
     */
    const double c0 = -1.0/30.0;
    const double c1 =  1.0/42.0;
    const double c2 = -1.0/30.0;
    const double ni2 = (1.0/n)*(1.0/n);
    const double ser =  ni2*ni2 * (c0 + ni2*(c1 + c2*ni2));
    result->val = (1.0 + 0.5/n + 1.0/(6.0*n*n) + ser) / n;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  return GSL_SUCCESS;
}


/* series for gammastar(x)
 * double-precision for x > 10.0
 */
static
int
gammastar_ser(const double x, gsl_sf_result * result)
{
  /* Use the Stirling series for the correction to Log(Gamma(x)),
   * which is better behaved and easier to compute than the
   * regular Stirling series for Gamma(x). 
   */
  const double y = 1.0/(x*x);
  const double c0 =  1.0/12.0;
  const double c1 = -1.0/360.0;
  const double c2 =  1.0/1260.0;
  const double c3 = -1.0/1680.0;
  const double c4 =  1.0/1188.0;
  const double c5 = -691.0/360360.0;
  const double c6 =  1.0/156.0;
  const double c7 = -3617.0/122400.0;
  const double ser = c0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7))))));
  result->val = exp(ser/x);
  result->err = 2.0 * GSL_DBL_EPSILON * result->val * GSL_MAX_DBL(1.0, ser/x);
  return GSL_SUCCESS;
}


/* x near a negative integer
 * Calculates sign as well as log(|gamma(x)|).
 * x = -N + eps
 * assumes N >= 1
 */
static
int
lngamma_sgn_sing(int N, double eps, gsl_sf_result * lng, double * sgn)
{
  if(eps == 0.0) {
    lng->val = 0.0;
    lng->err = 0.0;
    *sgn = 0.0;
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(N == 1) {
    /* calculate series for
     * g = eps gamma(-1+eps) + 1 + eps/2 (1+3eps)/(1-eps^2)
     * double-precision for |eps| < 0.02
     */
    const double c0 =  0.07721566490153286061;
    const double c1 =  0.08815966957356030521;
    const double c2 = -0.00436125434555340577;
    const double c3 =  0.01391065882004640689;
    const double c4 = -0.00409427227680839100;
    const double c5 =  0.00275661310191541584;
    const double c6 = -0.00124162645565305019;
    const double c7 =  0.00065267976121802783;
    const double c8 = -0.00032205261682710437;
    const double c9 =  0.00016229131039545456;
    const double g5 = c5 + eps*(c6 + eps*(c7 + eps*(c8 + eps*c9)));
    const double g  = eps*(c0 + eps*(c1 + eps*(c2 + eps*(c3 + eps*(c4 + eps*g5)))));

    /* calculate eps gamma(-1+eps), a negative quantity */
    const double gam_e = g - 1.0 - 0.5*eps*(1.0+3.0*eps)/(1.0 - eps*eps);

    lng->val = log(fabs(gam_e)/fabs(eps));
    lng->err = 2.0 * GSL_DBL_EPSILON * fabs(lng->val);
    *sgn = ( eps > 0.0 ? -1.0 : 1.0 );
    return GSL_SUCCESS;
  }
  else {
    double g;

    /* series for sin(Pi(N+1-eps))/(Pi eps) modulo the sign
     * double-precision for |eps| < 0.02
     */
    const double cs1 = -1.6449340668482264365;
    const double cs2 =  0.8117424252833536436;
    const double cs3 = -0.1907518241220842137;
    const double cs4 =  0.0261478478176548005;
    const double cs5 = -0.0023460810354558236;
    const double e2  = eps*eps;
    const double sin_ser = 1.0 + e2*(cs1+e2*(cs2+e2*(cs3+e2*(cs4+e2*cs5))));

    /* calculate series for ln(gamma(1+N-eps))
     * double-precision for |eps| < 0.02
     */
    double aeps = fabs(eps);
    double c1, c2, c3, c4, c5, c6, c7;
    double lng_ser;
    gsl_sf_result c0;
    gsl_sf_result psi_0;
    gsl_sf_result psi_1;
    gsl_sf_result psi_2;
    gsl_sf_result psi_3;
    gsl_sf_result psi_4;
    gsl_sf_result psi_5;
    gsl_sf_result psi_6;
    psi_2.val = 0.0;
    psi_3.val = 0.0;
    psi_4.val = 0.0;
    psi_5.val = 0.0;
    psi_6.val = 0.0;
    gsl_sf_lnfact_e(N, &c0);
    gsl_sf_psi_int_e(N+1, &psi_0);
    gsl_sf_psi_1_int_e(N+1, &psi_1);
    if(aeps > 0.00001) gsl_sf_psi_n_e(2, N+1.0, &psi_2);
    if(aeps > 0.0002)  gsl_sf_psi_n_e(3, N+1.0, &psi_3);
    if(aeps > 0.001)   gsl_sf_psi_n_e(4, N+1.0, &psi_4);
    if(aeps > 0.005)   gsl_sf_psi_n_e(5, N+1.0, &psi_5);
    if(aeps > 0.01)    gsl_sf_psi_n_e(6, N+1.0, &psi_6);
    c1 = psi_0.val;
    c2 = psi_1.val/2.0;
    c3 = psi_2.val/6.0;
    c4 = psi_3.val/24.0;
    c5 = psi_4.val/120.0;
    c6 = psi_5.val/720.0;
    c7 = psi_6.val/5040.0;
    lng_ser = c0.val-eps*(c1-eps*(c2-eps*(c3-eps*(c4-eps*(c5-eps*(c6-eps*c7))))));

    /* calculate
     * g = ln(|eps gamma(-N+eps)|)
     *   = -ln(gamma(1+N-eps)) + ln(|eps Pi/sin(Pi(N+1+eps))|)
     */
    g = -lng_ser - log(sin_ser);

    lng->val = g - log(fabs(eps));
    lng->err = c0.err + 2.0 * GSL_DBL_EPSILON * (fabs(g) + fabs(lng->val));

    *sgn = ( GSL_IS_ODD(N) ? -1.0 : 1.0 ) * ( eps > 0.0 ? 1.0 : -1.0 );

    return GSL_SUCCESS;
  }
  return GSL_SUCCESS;
}

/* x = eps near zero
 * gives double-precision for |eps| < 0.02
 */
static
int
lngamma_sgn_0(double eps, gsl_sf_result * lng, double * sgn)
{
  /* calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps/2 */
  const double c1  = -0.07721566490153286061;
  const double c2  = -0.01094400467202744461;
  const double c3  =  0.09252092391911371098;
  const double c4  = -0.01827191316559981266;
  const double c5  =  0.01800493109685479790;
  const double c6  = -0.00685088537872380685;
  const double c7  =  0.00399823955756846603;
  const double c8  = -0.00189430621687107802;
  const double c9  =  0.00097473237804513221;
  const double c10 = -0.00048434392722255893;
  const double g6  = c6+eps*(c7+eps*(c8 + eps*(c9 + eps*c10)));
  const double g   = eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*g6)))));

  /* calculate Gamma(eps) eps, a positive quantity */
  const double gee = g + 1.0/(1.0+eps) + 0.5*eps;

  lng->val = log(gee/fabs(eps));
  lng->err = 4.0 * GSL_DBL_EPSILON * fabs(lng->val);
  *sgn = GSL_SIGN(eps);

  return GSL_SUCCESS;
}


/* Lanczos method for real x > 0;
 * gamma=7, truncated at 1/(z+8) 
 * [J. SIAM Numer. Anal, Ser. B, 1 (1964) 86]
 */
static
int
lngamma_lanczos(double x, gsl_sf_result * result)
{
  int k;
  double Ag;
  double term1, term2;

  x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag = lanczos_7_c[0];
  for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

  /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
  term1 = (x+0.5)*log((x+7.5)/M_E);
  term2 = LogRootTwoPi_ + log(Ag);
  result->val  = term1 + (term2 - 7.0);
  result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(term1) + fabs(term2) + 7.0);
  result->err += GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}


inline
static
int
lngamma_1_pade(const double eps, gsl_sf_result * result)
{
  /* Use (2,2) Pade for Log[Gamma[1+eps]]/eps
   * plus a correction series.
   */
  const double n1 = -1.0017419282349508699871138440;
  const double n2 =  1.7364839209922879823280541733;
  const double d1 =  1.2433006018858751556055436011;
  const double d2 =  5.0456274100274010152489597514;
  const double num = (eps + n1) * (eps + n2);
  const double den = (eps + d1) * (eps + d2);
  const double pade = 2.0816265188662692474880210318 * num / den;
  const double c0 =  0.004785324257581753;
  const double c1 = -0.01192457083645441;
  const double c2 =  0.01931961413960498;
  const double c3 = -0.02594027398725020;
  const double c4 =  0.03141928755021455;
  const double eps5 = eps*eps*eps*eps*eps;
  const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
  result->val = eps * (pade + corr);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}


inline
static
int
lngamma_2_pade(const double eps, gsl_sf_result * result)
{
  /* Use (2,2) Pade for Log[Gamma[2+eps]]/eps
   * plus a correction series.
   */
  const double n1 = 1.000895834786669227164446568;
  const double n2 = 4.209376735287755081642901277;
  const double d1 = 2.618851904903217274682578255;
  const double d2 = 10.85766559900983515322922936;
  const double num = (eps + n1) * (eps + n2);
  const double den = (eps + d1) * (eps + d2);
  const double pade = 2.85337998765781918463568869 * num/den;
  const double c0 =  0.0001139406357036744;
  const double c1 = -0.0001365435269792533;
  const double c2 =  0.0001067287169183665;
  const double c3 = -0.0000693271800931282;
  const double c4 =  0.0000407220927867950;
  const double eps5 = eps*eps*eps*eps*eps;
  const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
  result->val = eps * (pade + corr);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}


/* gamma(x) for x >= 1/2
 * assumes x >= 1/2
 */
static
int
gamma_xgthalf(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.5) {
    result->val = 1.77245385090551602729817;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  } else if (x <= (GSL_SF_FACT_NMAX + 1.0) && x == floor(x)) {
    int n = (int) floor (x);
    result->val = fact_table[n - 1].f;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }    
  else if(fabs(x - 1.0) < 0.01) {
    /* Use series for Gamma[1+eps] - 1/(1+eps).
     */
    const double eps = x - 1.0;
    const double c1 =  0.4227843350984671394;
    const double c2 = -0.01094400467202744461;
    const double c3 =  0.09252092391911371098;
    const double c4 = -0.018271913165599812664;
    const double c5 =  0.018004931096854797895;
    const double c6 = -0.006850885378723806846;
    const double c7 =  0.003998239557568466030;
    result->val = 1.0/x + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*c7))))));
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(fabs(x - 2.0) < 0.01) {
    /* Use series for Gamma[1 + eps].
     */
    const double eps = x - 2.0;
    const double c1 =  0.4227843350984671394;
    const double c2 =  0.4118403304264396948;
    const double c3 =  0.08157691924708626638;
    const double c4 =  0.07424901075351389832;
    const double c5 = -0.00026698206874501476832;
    const double c6 =  0.011154045718130991049;
    const double c7 = -0.002852645821155340816;
    const double c8 =  0.0021039333406973880085;
    result->val = 1.0 + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 5.0) {
    /* Exponentiating the logarithm is fine, as
     * long as the exponential is not so large
     * that it greatly amplifies the error.
     */
    gsl_sf_result lg;
    lngamma_lanczos(x, &lg);
    result->val = exp(lg.val);
    result->err = result->val * (lg.err + 2.0 * GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }
  else if(x < 10.0) {
    /* This is a sticky area. The logarithm
     * is too large and the gammastar series
     * is not good.
     */
    const double gamma_8 = 5040.0;
    const double t = (2.0*x - 15.0)/5.0;
    gsl_sf_result c;
    cheb_eval_e(&gamma_5_10_cs, t, &c);
    result->val  = exp(c.val) * gamma_8;
    result->err  = result->val * c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x < GSL_SF_GAMMA_XMAX) {
    /* We do not want to exponentiate the logarithm
     * if x is large because of the inevitable
     * inflation of the error. So we carefully
     * use pow() and exp() with exact quantities.
     */
    double p = pow(x, 0.5*x);
    double e = exp(-x);
    double q = (p * e) * p;
    double pre = M_SQRT2 * M_SQRTPI * q/sqrt(x);
    gsl_sf_result gstar;
    int stat_gs = gammastar_ser(x, &gstar);
    result->val = pre * gstar.val;
    result->err = (x + 2.5) * GSL_DBL_EPSILON * result->val;
    return stat_gs;
  }
  else {
    OVERFLOW_ERROR(result);
  }
  return GSL_SUCCESS;
}


int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double * sgn)
{
  if(fabs(x - 1.0) < 0.01) {
    int stat = lngamma_1_pade(x - 1.0, result_lg);
    result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
    *sgn = 1.0;
    return stat;
  }
  else if(fabs(x - 2.0) < 0.01) {
   int stat = lngamma_2_pade(x - 2.0, result_lg);
    result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
    *sgn = 1.0;
    return stat;
  }
  else if(x >= 0.5) {
    *sgn = 1.0;
    return lngamma_lanczos(x, result_lg);
  }
  else if(x == 0.0) {
    *sgn = 0.0;
    DOMAIN_ERROR(result_lg);
  }
  else if(fabs(x) < 0.02) {
    return lngamma_sgn_0(x, result_lg, sgn);
  }
  else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
   /* Try to extract a fractional
     * part from x.
     */
    double z = 1.0 - x;
    double s = sin(M_PI*x);
    double as = fabs(s);
    if(s == 0.0) {
      *sgn = 0.0;
      DOMAIN_ERROR(result_lg);
    }
    else if(as < M_PI*0.015) {
      /* x is near a negative integer, -N */
      if(x < INT_MIN + 2.0) {
        result_lg->val = 0.0;
        result_lg->err = 0.0;
        *sgn = 0.0;
        GSL_ERROR ("error", GSL_EROUND);
      }
      else {
        int N = -(int)(x - 0.5);
        double eps = x + N;
        return lngamma_sgn_sing(N, eps, result_lg, sgn);
      }
    }
    else {
      gsl_sf_result lg_z;
      lngamma_lanczos(z, &lg_z);
      *sgn = (s > 0.0 ? 1.0 : -1.0);
      result_lg->val = M_LNPI - (log(as) + lg_z.val);
      result_lg->err = 2.0 * GSL_DBL_EPSILON * fabs(result_lg->val) + lg_z.err;
      return GSL_SUCCESS;
    }
  }
  else {
    /* |x| was too large to extract any fractional part */
    result_lg->val = 0.0;
    result_lg->err = 0.0;
    *sgn = 0.0;
    GSL_ERROR ("x too large to extract fraction part", GSL_EROUND);
  }
  return GSL_SUCCESS;
}


int gsl_sf_lngamma_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x - 1.0) < 0.01) {
    /* Note that we must amplify the errors
     * from the Pade evaluations because of
     * the way we must pass the argument, i.e.
     * writing (1-x) is a loss of precision
     * when x is near 1.
     */
    int stat = lngamma_1_pade(x - 1.0, result);
    result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
    return stat;
  }
  else if(fabs(x - 2.0) < 0.01) {
    int stat = lngamma_2_pade(x - 2.0, result);
    result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
    return stat;
  }
  else if(x >= 0.5) {
    return lngamma_lanczos(x, result);
  }
  else if(x == 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 0.02) {
    double sgn;
    return lngamma_sgn_0(x, result, &sgn);
  }
  else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
    /* Try to extract a fractional
     * part from x.
     */
    double z  = 1.0 - x;
    double s  = sin(M_PI*z);
    double as = fabs(s);
    if(s == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(as < M_PI*0.015) {
      /* x is near a negative integer, -N */
      if(x < INT_MIN + 2.0) {
        result->val = 0.0;
        result->err = 0.0;
        GSL_ERROR ("error", GSL_EROUND);
      }
      else {
        int N = -(int)(x - 0.5);
        double eps = x + N;
        double sgn;
        return lngamma_sgn_sing(N, eps, result, &sgn);
      }
    }
    else {
      gsl_sf_result lg_z;
      lngamma_lanczos(z, &lg_z);
      result->val = M_LNPI - (log(as) + lg_z.val);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg_z.err;
      return GSL_SUCCESS;
    }
  }
  else {
    /* |x| was too large to extract any fractional part */
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EROUND);
  }
  return GSL_SUCCESS;
}


int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= GSL_SF_FACT_NMAX){
    result->val = log(fact_table[n].f);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_lngamma_e(n+1.0, result);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result nf;
    gsl_sf_result mf;
    gsl_sf_result nmmf;
    if(m*2 > n) m = n-m;
    gsl_sf_lnfact_e(n, &nf);
    gsl_sf_lnfact_e(m, &mf);
    gsl_sf_lnfact_e(n-m, &nmmf);
    result->val  = nf.val - mf.val - nmmf.val;
    result->err  = nf.err + mf.err + nmmf.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  return GSL_SUCCESS;
}


int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if (n <= GSL_SF_FACT_NMAX) {
    result->val = (fact_table[n].f / fact_table[m].f) / fact_table[n-m].f;
    result->err = 6.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  } else {
    if(m*2 < n) m = n-m;

    if (n - m < 64)  /* compute product for a manageable number of terms */
      {
        double prod = 1.0;
        unsigned int k;
        
        for(k=n; k>=m+1; k--) {
          double tk = (double)k / (double)(k-m);
          if(tk > GSL_DBL_MAX/prod) {
            OVERFLOW_ERROR(result);
          }
          prod *= tk;
        }
        result->val = prod;
        result->err = 2.0 * GSL_DBL_EPSILON * prod * fabs(n-m);
        return GSL_SUCCESS;
      }
    else
      {
        gsl_sf_result lc;
        const int stat_lc = gsl_sf_lnchoose_e (n, m, &lc);
        const int stat_e  = gsl_sf_exp_err_e(lc.val, lc.err, result);
        return GSL_ERROR_SELECT_2(stat_lc, stat_e);
      }
  }
  return GSL_SUCCESS;
}


int
gsl_sf_gamma_e(const double x, gsl_sf_result * result)
{
  if(x < 0.5) {
    int rint_x = (int)floor(x+0.5);
    double f_x = x - rint_x;
    double sgn_gamma = ( GSL_IS_EVEN(rint_x) ? 1.0 : -1.0 );
    double sin_term = sgn_gamma * sin(M_PI * f_x) / M_PI;

    if(sin_term == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(x > -169.0) {
      gsl_sf_result g;
      gamma_xgthalf(1.0-x, &g);
      if(fabs(sin_term) * g.val * GSL_DBL_MIN < 1.0) {
        result->val  = 1.0/(sin_term * g.val);
        result->err  = fabs(g.err/g.val) * fabs(result->val);
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
      else {
        UNDERFLOW_ERROR(result);
      }
    }
    else {
      /* It is hard to control it here.
       * We can only exponentiate the
       * logarithm and eat the loss of
       * precision.
       */
      gsl_sf_result lng;
      double sgn;
      int stat_lng = gsl_sf_lngamma_sgn_e(x, &lng, &sgn);
      int stat_e   = gsl_sf_exp_mult_err_e(lng.val, lng.err, sgn, 0.0, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_lng);
    }
  }
  else {
    return gamma_xgthalf(x, result);
  }
  return GSL_SUCCESS;
}


int gsl_sf_fact_e (const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 18) {
    result->val = fact_table[n].f;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n <= GSL_SF_FACT_NMAX){
    result->val = fact_table[n].f;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
  return (1);
}


int gsl_sf_choose(unsigned int n, unsigned int m, double *result)
{
   gsl_sf_result gsl_result;
   int status = gsl_sf_choose_e(n, m, &gsl_result);
   if (status == GSL_SUCCESS) *result = gsl_result.val;
   return (status);
}


int gsl_sf_gamma(const double x, double *result)
{
   gsl_sf_result gsl_result;
   int status = gsl_sf_gamma_e (x, &gsl_result);
   if (status == GSL_SUCCESS) *result = gsl_result.val;
   return (status);
}

double gsl_sf_fact(const int n)
{
	if (n <= 0) return (1);
	else if(n <= GSL_SF_FACT_NMAX) {
		return (fact_table[n].f);
	} else {
		return (GSL_POSINF);
	}
}


