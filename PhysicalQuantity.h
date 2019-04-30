/*
 *  Created on: 2014-02-13
 *      Author: Jonas Östlund <uppfinnarjonas@gmail.com>
 *
 *  Extended on 2019-04-30 by Julien Pilet <julien.pilet@gmail.com>
 *
 *  HOW TO USE THESE CLASSES
 *  ========================
 *  These classes are useful in public APIs for
 *    * Parameter passing and return values of public functions
 *    * Instance variables of objects
 *
 *  In the hidden or low-level parts of the implementation,
 *  it may be more convenient to use raw floating point values.
 */


#ifndef PHYSICALQUANTITY_H_
#define PHYSICALQUANTITY_H_

#include <type_traits>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

namespace physical_quantity {

#define QUANTITY_TEMPLATE \
	 int TimeDim/*t*/, int LengthDim/*l*/, int AngleDim/*a*/, int MassDim/*m*/, \
	 int TemperatureDim /*T*/, int ElectricCurrentDim /*I*/, int AmountOfSubstanceDim /*N*/, int LuminousIntensityDim /*J*/

#define QUANTITY_ARGS TimeDim, LengthDim, AngleDim, MassDim, TemperatureDim, ElectricCurrentDim, AmountOfSubstanceDim, LuminousIntensityDim

template <typename T, typename System, QUANTITY_TEMPLATE> class PhysicalQuantity;


#define FOREACH_SI_PREFIX(OP, type, name) \
	OP(type, yotta##name, 1e24, 0) \
	OP(type, zetta##name, 1e21, 0) \
	OP(type, exa##name, 1e18, 0) \
	OP(type, peta##name, 1e15, 0) \
	OP(type, tera##name, 1e12, 0) \
	OP(type, giga##name, 1e9, 0) \
	OP(type, mega##name, 1e6, 0) \
	OP(type, kilo##name, 1e3, 0) \
	OP(type, name, 1.0, 0) \
	OP(type, milli##name, 1e-3, 0) \
	OP(type, micro##name, 1e-6, 0) \
	OP(type, nano##name, 1e-9, 0) \
	OP(type, pico##name, 1e-12, 0) \
	OP(type, femto##name, 1e-15, 0) \
	OP(type, atto##name, 1e-18, 0) \
	OP(type, zepto##name, 1e-21, 0) \
	OP(type, yocto##name, 1e-24, 0)

#define PREFIX_TIME(name, factor) OP(Time, name ## seconds, factor)

// OP(type, name, factor)
#define FOREACH_TIME_UNIT(OP) \
    FOREACH_SI_PREFIX(OP, Time, seconds) \
    OP(Time, minutes, 60.0, 0) \
    OP(Time, hours, 3600.0, 0) \
    OP(Time, days, 24*3600.0, 0) \
    OP(Time, weeks, 7*24*3600.0, 0)

#define FOREACH_LENGTH_UNIT(OP) \
  FOREACH_SI_PREFIX(OP, Length, meters) \
  OP(Length, nauticalMiles, 1852.0, 0) // 1 nautical mile = 1852.0 m

#define FOREACH_ANGLE_UNIT(OP) \
  OP(Angle, radians, 1.0, 0) \
  OP(Angle, degrees, 6.28318530718/360.0, 0) // 1 degree = M_PI/180 radians <=> 180 degrees = M_PI radians

#define FOREACH_VELOCITY_UNIT(OP) \
  OP(Velocity, metersPerSecond, 1.0, 0) \
  OP(Velocity, knots, 1852.0/3600.0, 0) \
  OP(Velocity, kilometersPerHour, 1000.0/3600.0, 0) \
  OP(Velocity, milesPerHour, 1609.0/3600.0, 0)

#define FOREACH_ACCELERATION_UNIT(OP) \
  OP(Acceleration, metersPerSecondSquared, 1.0, 0)

#define FOREACH_MASS_UNIT(OP) \
  FOREACH_SI_PREFIX(OP, Mass, grams) \
  OP(Mass, skeppund, 170.0 * 1000, 0) \
  OP(Mass, lispund, (170.0 * 1000)/20, 0)

#define FOREACH_ANGULAR_VELOCITY_UNIT(OP) \
  OP(AngularVelocity, radiansPerSecond, 1.0, 0) \
  OP(AngularVelocity, degreesPerSecond, 6.28318530718/360.0, 0) \
  OP(AngularVelocity, rpm, 6.28318530718 / 60.0, 0) // 1 rpm = 2 * pi rad in 60 sec.

#define FOREACH_ELECTRIC_CURRENT_UNIT(OP) \
  FOREACH_SI_PREFIX(OP, ElectricCurrent, ampere)

#define FOREACH_TEMPERATURE_UNIT(OP) \
  OP(Temperature, kelvin, 1, 0) \
  OP(Temperature, celsius, 1, 273.15)

#define FOREACH_AMOUNT_OF_SUBSTANCE_UNIT(OP) \
  OP(AmountOfSubstance, mole, 1, 0)

#define FOREACH_LUMINOUS_INTENSITY_UNIT(OP) \
  OP(LuminousIntensity, candela, 1, 0)

#define FOREACH_UNIT(OP) \
  FOREACH_TIME_UNIT(OP) \
  FOREACH_LENGTH_UNIT(OP) \
  FOREACH_ANGLE_UNIT(OP) \
  FOREACH_VELOCITY_UNIT(OP) \
  FOREACH_ACCELERATION_UNIT(OP) \
  FOREACH_MASS_UNIT(OP) \
  FOREACH_ANGULAR_VELOCITY_UNIT(OP) \
  FOREACH_ELECTRIC_CURRENT_UNIT(OP) \
  FOREACH_TEMPERATURE_UNIT(OP) \
  FOREACH_AMOUNT_OF_SUBSTANCE_UNIT(OP) \
  FOREACH_LUMINOUS_INTENSITY_UNIT(OP)

#define FOREACH_QUANTITY(OP) \
  OP(Time, seconds,              1, 0, 0, 0, 0, 0, 0, 0) \
  OP(Length, meters,             0, 1, 0, 0, 0, 0, 0, 0) \
  OP(Angle, radians,             0, 0, 1, 0, 0, 0, 0, 0) \
  OP(Mass, kilograms,            0, 0, 0, 1, 0, 0, 0, 0) \
  OP(Temperature, kelvin,        0, 0, 0, 0, 1, 0, 0, 0) \
  OP(ElectricCurrent, ampere,    0, 0, 0, 0, 0, 1, 0, 0) \
  OP(AmountOfSubstance, mole,    0, 0, 0, 0, 0, 0, 1, 0) \
  OP(LuminousIntensity, candela, 0, 0, 0, 0, 0, 0, 0, 1) \
  OP(AngularVelocity, radiansPerSecond, -1, 0, 1, 0, 0, 0, 0, 0) \
  OP(Velocity, metersPerSecond, -1, 1, 0, 0, 0, 0, 0, 0) \
  OP(Acceleration, metersPerSecondSquared, -2, 1, 0, 0, 0, 0, 0, 0)

enum class Quantity {
  // Any quantity that has not been declared maps to this one.
  Undeclared,
#define MAKE_QUANTITY_ENUM(name, siUnit, t, l, a, m, T, I, N, J) \
    name,
FOREACH_QUANTITY(MAKE_QUANTITY_ENUM)
#undef MAKE_QUANTITY_ENUM
};

enum class Unit {
  // This unit can be used with any quantity and is
  // the SI unit of that quantity. So even if we are working
  // with volume and there are no volume quantities declared,
  // this unit will in that case represents cubic meter.
  SIUnit,
#define MAKE_UNIT_ENUM(type, name, factor, offset) \
  name,
FOREACH_UNIT(MAKE_UNIT_ENUM)
#undef MAKE_UNIT_ENUM
};

/*
 * This general template will, in particular, apply to
 * to Unit::SIUnit, which is the default unit for all
 * quantities that we have not declared, such as acceleration.
 */
template <Unit i>
struct UnitInfo {
  static constexpr double getFactor() {return 1.0;}
  static constexpr double getOffset() {return 0.0;}
  static const Quantity quantity = Quantity::Undeclared;
  static const Unit unit = i;
};

// Here, the above template is being specialized for all declared units.
#define MAKE_INT_TO_UNIT_FACTOR(type, name, factor, offset) \
  template<> struct UnitInfo<Unit::name> { \
  static constexpr double getFactor() { return factor; } \
  static constexpr double getOffset() { return offset; } \
  static const Quantity quantity = Quantity::type; \
  static const Unit unit = Unit::name; \
};
FOREACH_UNIT(MAKE_INT_TO_UNIT_FACTOR)
#undef MAKE_INT_TO_UNIT_FACTOR

/*
 * Default template that catches all quantities
 * that we have not declared, but that are nevertheless
 * valid, such as acceleration.
 */
template <typename unitsys, QUANTITY_TEMPLATE>
struct QuantityInfo {
  static const Quantity quantity = Quantity::Undeclared;
  static const Unit unit = Unit::SIUnit;
};

#define MAKE_QUANTITY_INFO(name, siUnit, t, l, a, m, T, I, N, J) \
  template <typename unitsys> \
  struct QuantityInfo<unitsys, t, l, a, m, T, I, N, J> { \
    static const Unit unit = unitsys::name; \
    static const Quantity quantity = Quantity::name; \
  };
FOREACH_QUANTITY(MAKE_QUANTITY_INFO)
#undef MAKE_QUANTITY_INFO


// Converts between arbitrary units. This code is used
// by PhysicalQuantity.
template <typename T, Unit FromUnit, Unit ToUnit>
struct ConvertUnit {
  static_assert(
      UnitInfo<FromUnit>::quantity == UnitInfo<ToUnit>::quantity
      || FromUnit == Unit::SIUnit || ToUnit == Unit::SIUnit,
      "Incompatible quantitites");
  static T apply(T x) {
    constexpr double f = UnitInfo<FromUnit>::getFactor()/UnitInfo<ToUnit>::getFactor();
    constexpr double o = (UnitInfo<FromUnit>::getOffset() - UnitInfo<ToUnit>::getOffset()) / UnitInfo<ToUnit>::getFactor();
    return T(f)*x + T(o);
  }
};
// For the special case when the unit is the same,
// there is no need to multiply.
template <typename T, Unit u>
struct ConvertUnit<T, u, u> {
public:
  static T apply(T x) {return x;}
};

namespace UnitSystem {
  struct SI {
#define MAKE_SI_UNIT(name, siUnit, t, l, a, m, T, I, N, J) \
  static const Unit name = Unit::siUnit;
FOREACH_QUANTITY(MAKE_SI_UNIT)
#undef MAKE_SI_UNIT
  };

  // Suitable for integer and fixed-point representation
  struct CustomAnemoUnits : public SI {
    static const Unit Time = Unit::milliseconds;
    static const Unit Angle = Unit::degrees;
    static const Unit Velocity = Unit::knots;
  };
};

template <typename T, typename System, QUANTITY_TEMPLATE>
std::string physQuantToString(const PhysicalQuantity<T, System, QUANTITY_ARGS> &x);

struct NonsenseType {};

template <typename T, typename System, QUANTITY_TEMPLATE>
struct DimensionlessTraits {
  typedef NonsenseType PublicType;
  typedef T PrivateType;
  static NonsenseType get(T x) {return NonsenseType{};}
};

template <typename T, typename System>
struct DimensionlessTraits<T, System, 0, 0, 0, 0, 0, 0, 0, 0> {
  typedef NonsenseType PrivateType;
  typedef T PublicType;
  static T get(T x) {return x;}
};


template <typename T, typename System, QUANTITY_TEMPLATE>
class PhysicalQuantity {
public:
  typedef DimensionlessTraits<T, System, QUANTITY_ARGS> DimensionlessInfo;

  typedef T ValueType;
  typedef System SystemType;

  PhysicalQuantity() : _x(T(std::numeric_limits<double>::signaling_NaN())) {} // TODO: FIX THIS!!!

  typedef QuantityInfo<System, QUANTITY_ARGS> QInfo;
  typedef UnitInfo<QInfo::unit> UInfo;

  typedef PhysicalQuantity<T, System, QUANTITY_ARGS> ThisType;

  static constexpr bool isDimensionless = TimeDim == 0 && LengthDim == 0
      && AngleDim == 0 && MassDim == 0 && TemperatureDim == 0
      && ElectricCurrentDim == 0 && AmountOfSubstanceDim == 0 && LuminousIntensityDim == 0;

  // For example, for "Angle<T>" and "degrees", the following will create within
  // the class Angle<T>:
  // static Angle<T> degrees(T x);
  // static Angle<T> make_degrees(T x);
  // T degrees() const;
  //
  // make_degrees and degrees are synonms, but they help removing the ambiguity
  // when getting a function pointer: &Angle<T>::make_degrees has no ambiguity.
#define MAKE_UNIT_CONVERTERS(type, name, factor, offset) \
  static ThisType name(T x) { return make_##name(x); }\
  static ThisType make_##name(T x) { \
    static_assert(UnitInfo<Unit::name>::quantity == QInfo::quantity, "Incompatible unit and quantity"); \
    return ThisType(ConvertUnit<T, Unit::name, System::type>::apply(x)); \
  } \
  T name() const { \
    static_assert(UnitInfo<Unit::name>::quantity == QInfo::quantity, "Incompatible unit and quantity"); \
    return ConvertUnit<T, System::type, Unit::name>::apply(_x); \
  }
  FOREACH_UNIT(MAKE_UNIT_CONVERTERS)

  template <Unit unit>
  static ThisType make(T x) {
    static_assert(UnitInfo<unit>::quantity == QInfo::quantity
        || unit == Unit::SIUnit,
        "Incompatible unit and quantity");
    return ThisType(ConvertUnit<T, unit, UInfo::unit>::apply(x));
  }

  static ThisType makeFromSI(T x) {
    return make<Unit::SIUnit>(x);
  }

  template <Unit unit>
  T get() const {
    static_assert(UnitInfo<unit>::quantity == QInfo::quantity
        || unit == Unit::SIUnit,
        "Incompatible unit and quantity");
    return ConvertUnit<T, UInfo::unit, unit>::apply(_x);
  }

  T getSI() const {
    return get<Unit::SIUnit>();
  }

  static PhysicalQuantity<T, System, QUANTITY_ARGS> wrap(T x) {
    return PhysicalQuantity<T, System, QUANTITY_ARGS>(x);
  }

  ThisType operator*(T s) const {
    return ThisType(s*_x);
  }

  ThisType scaled(T s) const {
    return ThisType(s*_x);
  }

  static ThisType zero() {
    return ThisType(T(0.0));
  }

  T dimensionless() const {
    static_assert(isDimensionless, "Only applicable to Dimensionlesss");
    return ConvertUnit<T, UInfo::unit, Unit::SIUnit>::apply(_x);
  }

  static ThisType dimensionless(T x) {
    static_assert(isDimensionless, "Only applicable to Dimensionless types");
    return ThisType(ConvertUnit<T, Unit::SIUnit, UInfo::unit>::apply(x));
  }

  ThisType operator+(const ThisType &other) const {
    return ThisType(_x + other._x);
  }

  ThisType operator-(const ThisType &other) const {
    return ThisType(_x - other._x);
  }

  ThisType &operator+=(const ThisType &other) {
    _x += other._x;
    return *this;
  }

  ThisType &operator-=(const ThisType &other) {
    _x -= other._x;
    return *this;
  }

  ThisType fabs() const { return ThisType(::fabs(_x)); }

  bool isNaNQuantity() const { return physical_quantity::isNaN(_x); }

  bool isFiniteQuantity() const { return physical_quantity::isFinite(_x); }

  ThisType operator-() const {
    return ThisType(-_x);
  }

  bool operator < (ThisType other) const {return _x < other._x;}
  bool operator <= (ThisType other) const {return _x <= other._x;}
  bool operator > (ThisType other) const {return _x > other._x;}
  bool operator >= (ThisType other) const {return _x >= other._x;}
  bool operator == (ThisType other) const {return _x == other._x;}

  template <typename S>
  PhysicalQuantity<S, System, QUANTITY_ARGS> cast() const {
    typedef PhysicalQuantity<S, System, QUANTITY_ARGS> DstType;
    static constexpr Unit unit = DstType::UInfo::unit;
    return DstType::template make<unit>(static_cast<S>(get<unit>()));
  }

  template <typename S> // TODO: Maybe allow to cast between different systems.
  operator PhysicalQuantity<S, System, QUANTITY_ARGS>() const {
    return cast<S>();
  }

  void sincos(T *sinAngle, T *cosAngle) const {
    static_assert(UInfo::quantity == Quantity::Angle, "Only applicable to angles");
    T rad = radians();
    *sinAngle = sin(rad);
    *cosAngle = cos(rad);
  }

  std::string str() const {
    return physQuantToString(*this);
  }

  // Special method returning true for the comparison nan == nan.
  bool eqWithNan(ThisType other) const {
    return strictEquality(_x, other._x);
  }

  bool nearWithNan(ThisType other, double marg) const {
    return physical_quantity::nearWithNan(_x, other._x, marg);
  }

  operator typename DimensionlessInfo::PublicType() const {
    return DimensionlessInfo::get(_x);
  }
private:
  PhysicalQuantity(T x) : _x(x) {}
  T _x;
};

template <typename T, typename sys,
  int t0, int l0, int a0, int m0, int T0, int I0, int N0, int J0,
  int t1, int l1, int a1, int m1, int T1, int I1, int N1, int J1>
struct Division {
  typedef PhysicalQuantity<T, sys, t0 - t1, l0 - l1, a0 - a1, m0 - m1, T0 - T1, I0 - I1, N0 - N1, J0 - J1> DstType;
  typedef PhysicalQuantity<T, sys, t0, l0, a0, m0, T0, I0, N0, J0> NumeratorType;
  typedef PhysicalQuantity<T, sys, t1, l1, a1, m1, T1, I1, N1, J1> DenominatorType;

  static DstType apply(const NumeratorType &a, const DenominatorType &b) {
    T aValue = a.getSI();
    T bValue = b.getSI();
    return DstType::makeFromSI(aValue/bValue);
  }
};

template <typename Storage, typename sys,
  int t, int l, int a, int m, int T, int I, int N, int J>
struct Division<Storage, sys,
	t, l, a, m, T, I, N, J,
	t, l, a, m, T, I, N, J> {
  static PhysicalQuantity<Storage, sys, 0, 0, 0, 0, 0, 0, 0, 0>
    apply(const PhysicalQuantity<Storage, sys, t, l, a, m, T, I, N, J> &A,
          const PhysicalQuantity<Storage, sys, t, l, a, m, T, I, N, J> &B) {
    /*
     *  We have some code out there that divides
     *  two quantities like this, where T is fixed-point
     *  or integer. So this specialization avoids the unnecessary
     *  loss of precision due to first converting it to an SI unit.
     */
    static const Unit unit = PhysicalQuantity<Storage, sys, t, l, a, m>::UInfo::unit;
    T aValue = A.template get<unit>();
    T bValue = B.template get<unit>();
    typedef PhysicalQuantity<Storage, sys, 0, 0, 0, 0> DstType;
    return DstType::dimensionless(aValue/bValue);
  }
};

template <typename T, typename sys,
  int t0, int l0, int a0, int m0, int T0, int I0, int N0, int J0,
  int t1, int l1, int a1, int m1, int T1, int I1, int N1, int J1>
PhysicalQuantity<T, sys, t0 - t1, l0 - l1, a0 - a1, m0 - m1, T0 - T1, I0 - I1, N0 - N1, J0 - J1>
operator/(
    const PhysicalQuantity<T, sys, t0, l0, a0, m0, T0, I0, N0, J0> &a,
    const PhysicalQuantity<T, sys, t1, l1, a1, m1, T1, I1, N1, J1> &b) {
  return Division<T, sys, t0, l0, a0, m0, T0, I0, N0, J0, t1, l1, a1, m1, T1, I1, N1, J1>::apply(a, b);
}

template <typename T, typename sys,
  int t0, int l0, int a0, int m0, int T0, int I0, int N0, int J0,
  int t1, int l1, int a1, int m1, int T1, int I1, int N1, int J1>
static PhysicalQuantity<T, sys, t0 + t1, l0 + l1, a0 + a1, m0 + m1, T0 + T1, I0 + I1, N0 + N1, J0 + J1>
  operator*(
    const PhysicalQuantity<T, sys, t0, l0, a0, m0, T0, I0, N0, J0> &a,
    const PhysicalQuantity<T, sys, t1, l1, a1, m1, T1, I1, N1, J1> &b) {
  T aValue = a.getSI();
  T bValue = b.getSI();
  typedef PhysicalQuantity<T, sys, t0 + t1, l0 + l1, a0 + a1, m0 + m1, T0 + T1, I0 + I1, N0 + N1, J0 + J1> DstType;
  return DstType::makeFromSI(aValue*bValue);
}

template <typename A, typename B>
using Per = decltype((std::declval<A>()) / (std::declval<B>()));

template <typename T>
using TimeDerivative = Per<T,
    PhysicalQuantity<typename T::ValueType, typename T::SystemType, 1, 0, 0, 0, 0, 0, 0, 0>>;

// http://en.cppreference.com/w/cpp/language/type_alias
template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Dimensionless = PhysicalQuantity<T, System, 0, 0, 0, 0, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Duration = PhysicalQuantity<T, System, 1, 0, 0, 0, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Frequency = PhysicalQuantity<T, System, -1, 0, 0, 0, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Length = PhysicalQuantity<T, System, 0, 1, 0, 0, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Velocity = PhysicalQuantity<T, System, -1, 1, 0, 0, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Acceleration = PhysicalQuantity<T, System, -2, 1, 0, 0, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Angle = PhysicalQuantity<T, System, 0, 0, 1, 0, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Mass = PhysicalQuantity<T, System, 0, 0, 0, 1, 0, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using Temperature = PhysicalQuantity<T, System, 0, 0, 0, 0, 1, 0, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using ElectricCurrent = PhysicalQuantity<T, System, 0, 0, 0, 0, 0, 1, 0, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using AmountOfSubstance = PhysicalQuantity<T, System, 0, 0, 0, 0, 0, 0, 1, 0>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using LuminousIntensity = PhysicalQuantity<T, System, 0, 0, 0, 0, 0, 0, 0, 1>;

template <typename T=double, typename System=UnitSystem::CustomAnemoUnits>
using AngularVelocity = TimeDerivative<Angle<T, System>>;

template <typename T, typename System, QUANTITY_TEMPLATE>
PhysicalQuantity<T, System, QUANTITY_ARGS> operator*(T s,
    const PhysicalQuantity<T, System, QUANTITY_ARGS> &x) {
  return x*s;
}

template <typename T, typename sys, QUANTITY_TEMPLATE>
std::string physQuantToString(const PhysicalQuantity<T, sys, QUANTITY_ARGS> &x) {
  std::stringstream ss;

  const auto exponent = [](int e) -> std::string {
    if (e) {
      return std::string("^") + std::to_string(e);
    }
    return "";
  };

  ss << x.getSI() << " [";
  if (TimeDim != 0) { ss << " s" << exponent(TimeDim); }
  if (LengthDim != 0) { ss << " m" << exponent(LengthDim); }
  if (AngleDim != 0) { ss << " r" << exponent(AngleDim); }
  if (MassDim != 0) { ss << " g" << exponent(MassDim); }
  if (TemperatureDim) { ss << " K" << exponent(TemperatureDim); }
  if (ElectricCurrentDim) { ss << " A" << exponent(ElectricCurrentDim); }
  if (AmountOfSubstanceDim) { ss << " N" << exponent(AmountOfSubstanceDim); }
  if (LuminousIntensityDim) { ss << " J" << exponent(LuminousIntensityDim); }

  ss << " ]";
  return ss.str();
}

template <typename T, typename System>
std::string physQuantToString(const Duration<T, System> &x) {
   std::stringstream ss;
   Duration<T> remaining(x);
#define FORMAT_DURATION_UNIT(unit) \
   if (remaining.unit() >= 1) { \
     if (ss.str().size() > 0) {ss << ", ";} \
     ss << floor(remaining.unit()) << " " #unit ; \
     remaining -= Duration<T>::unit(floor(remaining.unit())); \
   }
   FORMAT_DURATION_UNIT(weeks)
   FORMAT_DURATION_UNIT(days)
   FORMAT_DURATION_UNIT(hours)
   FORMAT_DURATION_UNIT(minutes)
   FORMAT_DURATION_UNIT(seconds)
#undef FORMAT_DURATION_UNIT
   return ss.str();
}

/////////////////////////////////////////////////////// Extra stuff
template<typename T, int N>
class FixedArray {
 public:
  FixedArray() { }
  FixedArray(const FixedArray& a) {
    for (int i = 0; i < N; ++i) {
      _data[i] = a._data[i];
    }
  }

  T& operator[](int i) { return _data[i]; }
  const T& operator[](int i) const { return _data[i]; }

  const T* data() const { return _data; }
 private:
  T _data[N];
};

template <typename T, int N>
class Vectorize : public FixedArray<T, N> {
  public:
    Vectorize<T, N>(std::initializer_list<T> list) {
        int i=0;
        for (T element : list) {
            (*this)[i++] = element;
        }
    }

    explicit Vectorize<T, N>(const T x[N]) {
      for (int i = 0; i < N; i++) {
        (*this)[i] = x[i];
      }
    }

    static Vectorize<T, N> all(T value) {
        Vectorize<T, N> result;
        for (int i = 0; i < N; ++i) {
            result[i] = value;
        }
        return result;
    }

    Vectorize<T, N> operator + (const Vectorize<T, N>& other) const {
        Vectorize result;
        for (int i = 0; i < N; ++i) {
            result[i] = (*this)[i] + other[i];
        }
        return result;
    }

    Vectorize<T, N> operator - (const Vectorize<T, N>& other) const {
        Vectorize result;
        for (int i = 0; i < N; ++i) {
            result[i] = (*this)[i] - other[i];
        }
        return result;
    }

    template <typename FactorType>
    Vectorize<T, N> scaled(FactorType factor) const {
        Vectorize result;
        for (int i = 0; i < N; ++i) {
            result[i] = (*this)[i].scaled(factor);
        }
        return result;
    }

    Vectorize<T, N> operator- () const {
      return Vectorize<T, N>{-(*this)[0], -(*this)[1]};
    }

    bool operator == (const Vectorize<T, N>& other) const {
        for (int i = 0; i < N; ++i) {
          if (!((*this)[i] == other[i])) return false;
        }
        return true;
    }

    Vectorize<T, N>& operator+= (const Vectorize<T, N>& other) {
      for (int i = 0; i < N; ++i) {
        (*this)[i] += other[i];
      }
      return *this;
    }

    // General purpose implementation (doesn't work with FixedPoint I think)
    T norm() const {
      constexpr Unit unit = T::UInfo::unit;
      typedef typename T::ValueType type;
      type sum(0.0);
      for (int i = 0; i < N; i++) {
        type x = ((*this)[i]).template get<unit>();
        sum += x*x;
      }
      return T::template make<unit>(sqrt(sum));
    }

    Vectorize() { }
  private:
};

template <typename FactorType, typename ElemType, int N>
Vectorize<ElemType, N> operator*(FactorType x, const Vectorize<ElemType, N> &v) {
  return v.scaled(x);
}

}  // namespace physical_quantity

// Overloading sin and cos must be done outside namespace
template <typename T, typename System>
T cos(physical_quantity::Angle<T, System> x) {return cos(x.radians());}

template <typename T, typename System>
T sin(physical_quantity::Angle<T, System> x) {return sin(x.radians());}

template <typename T, typename System>
T tan(physical_quantity::Angle<T, System> x) {return tan(x.radians());}

template <typename T, typename s, QUANTITY_TEMPLATE>
physical_quantity::PhysicalQuantity<T, s, QUANTITY_ARGS> fabs(physical_quantity::PhysicalQuantity<T, s, QUANTITY_ARGS> x) {
  return x.fabs();
}

namespace physical_quantity {

template <typename T, typename System = UnitSystem::CustomAnemoUnits>
class HorizontalMotion : public Vectorize<Velocity<T, System>, 2> {
  public:
    typedef Velocity<T, System> InnerType;
    typedef Vectorize<Velocity<T, System>, 2> BaseType;

    HorizontalMotion(InnerType eastWest, InnerType southNorth) {
      (*this)[0] = eastWest;
      (*this)[1] = southNorth;
    }

    HorizontalMotion(const BaseType& base) : BaseType(base) { }

    HorizontalMotion() { }

    static HorizontalMotion<T, System> zero() {
        return HorizontalMotion(BaseType::all(InnerType::zero()));
    }

    static HorizontalMotion<T, System> polar(Velocity<T, System> speed, Angle<T, System> direction) {
        // A direction of 0 points to north.
        T sinDir, cosDir;
        direction.sincos(&sinDir, &cosDir);
        return HorizontalMotion<T>(
            speed.scaled(sinDir),
            speed.scaled(cosDir));
    }

    // Special implementation, that also works with fixed point, I think.
    Velocity<T, System> norm() const {
        T a = (*this)[0].knots();
        T b = (*this)[1].knots();
        return Velocity<T>::knots(sqrt(a*a + b*b));
    }

    Angle<T, System> angle() const {
        return Angle<T, System>::radians(atan2(
                (*this)[0].knots(),
                (*this)[1].knots()));
    }

    template <typename Dst>
    HorizontalMotion<Dst, System> cast() const {
      return HorizontalMotion<Dst, System>((*this)[0], (*this)[1]);
    }

    // Define what the vector dimensions mean.
    enum {
        EAST_TO_WEST = 0,
        SOUTH_TO_NORTH = 1
    };
};

template <typename T, typename s, QUANTITY_TEMPLATE>
bool isFinite(const PhysicalQuantity<T, s, QUANTITY_ARGS> &x) {
  return x.isFiniteQuantity();
}
template <typename T, typename s, QUANTITY_TEMPLATE>
bool isNaN(const PhysicalQuantity<T, s, QUANTITY_ARGS> &x) {
  return x.isNaNQuantity();
}

template <typename T>
// If T is an integral type (char, int long, etc), the following type resolves
// as bool and is used as return type. If the return type can't be deduced from
// T, this function will not be used for this type (SFINAE)
typename std::enable_if<std::is_integral<T>::value, bool>::type isNaN(T x) {
	return false;
}

template <typename T>
typename std::enable_if<std::is_integral<T>::value, bool>::type isFinite(T x) {
	return true;
}

template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type isNaN(T x) {
	return std::isnan(x);
}

template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type isFinite(T x) {
	return std::isfinite(x);
}


// Literals
#define DEFINE_LITERAL(QUANTITY, WHAT, LIT) \
  inline QUANTITY<double> operator"" LIT (long double x) { \
    return QUANTITY<double>::WHAT(x); \
}

DEFINE_LITERAL(Angle, degrees, _deg)
DEFINE_LITERAL(Angle, radians, _rad)
DEFINE_LITERAL(Length, meters, _m)
DEFINE_LITERAL(Length, kilometers, _km)
DEFINE_LITERAL(Length, nauticalMiles, _M)
DEFINE_LITERAL(Duration, seconds, _s)
DEFINE_LITERAL(Duration, seconds, _seconds)
DEFINE_LITERAL(Duration, days, _days)
DEFINE_LITERAL(Duration, hours, _h)
DEFINE_LITERAL(Duration, hours, _hours)
DEFINE_LITERAL(Duration, minutes, _minutes)
DEFINE_LITERAL(Velocity, metersPerSecond, _mps)
DEFINE_LITERAL(Velocity, knots, _kn)
DEFINE_LITERAL(Velocity, knots, _kt)
DEFINE_LITERAL(Acceleration, metersPerSecondSquared, _mps2)
#undef DEFINE_LITERAL

}  // namespace physical_quantity

#endif /* PHYSICALQUANTITY_H_ */
