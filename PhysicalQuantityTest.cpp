/*
 *  Created on: 2014-02-13
 *      Author: Jonas Östlund <uppfinnarjonas@gmail.com>
 */

#include <cmath>
#include "PhysicalQuantity.h"
#include <gtest/gtest.h>

using namespace physical_quantity;


TEST(PhysQuantTest, ScaleByDimensionless) {
  Dimensionless<double> factor = Dimensionless<double>::dimensionless(2.0);
  Angle<double> angle = Angle<double>::degree(34.5);
  Angle<double> product = factor*angle;
  EXPECT_NEAR(product.degree(), 69.0, 1.0e-6);

  Angle<double> divided = product/factor;
  EXPECT_NEAR(divided.degree(), 34.5, 1.0e-6);
}

TEST(PhysQuantTest, DimensionlessTest) {
  Dimensionless<double> x = Dimensionless<double>::dimensionless(34.4);
  double y = x;
  EXPECT_EQ(y, 34.4);
}

TEST(PhysQuantTest, MixingQuantities) {
  Velocity<double> vel = Velocity<double>::knot(2.4);
  EXPECT_NEAR(vel.meterPerSecond(), 1.2346666666666667, 1.0e-6);
  // 2.4 knot = 1.23466667 meter per second

  Duration<double> dur = Duration<double>::hour(3.4);
  EXPECT_NEAR(dur.second(), 12240, 1.0e-6);
  // 3.4 hours = 12240 second

  Length<double> dist = Length<double>::nauticalMile(8.16);
  EXPECT_NEAR(dist.meter(), 15112.32, 1.0e-6);
  // 15112.32 meter

  Length<double> computedDist = vel*dur;
  EXPECT_NEAR(computedDist.meter(), 15112.3200408, 0.1);

  Velocity<double> computedVelocity = dist/dur;
  EXPECT_NEAR(computedVelocity.meterPerSecond(), 1.23466667, 0.001);

  Duration<double> computedDuration = dist/vel;
  EXPECT_NEAR(computedDuration.second(), 12240, 0.1);
}

TEST(PhysQuantTest, CircumferenceTest) {
  Length<double> circumference = Length<double>::meter(4.0e7);
  double minutes = 60.0*Angle<double>::radian(6.28318530718).degree();
  double oneNauticalMileMeters = circumference.meter()/minutes;
  EXPECT_NEAR(oneNauticalMileMeters, Length<double>::nauticalMile(1.0).meter(), 30.0);
  EXPECT_NEAR(circumference.nauticalMile(), minutes, 40);
}

TEST(PhysQuantTest, DurationTest) {
  double n = 34.0;
  EXPECT_NEAR(Duration<double>::second(n).second(), n, 1.0e-9);
  EXPECT_NEAR(Duration<double>::second(60).minute(), 1.0, 1.0e-9);
}

TEST(PhysQuantTest, HydroptereTest) {
  // « sustaining a speed of 52.86 knot (97.90 km/h; 60.83 mph) »
  Velocity<double> va = Velocity<double>::knot(52.86);
  Velocity<double> vb = Velocity<double>::kilometerPerHour(97.90);
  Velocity<double> vc = Velocity<double>::milePerHour(60.83);
  EXPECT_NEAR(va.meterPerSecond(), vb.meterPerSecond(), 0.1);
  EXPECT_NEAR(va.meterPerSecond(), vc.meterPerSecond(), 0.1);
}

TEST(PhysQuantTest, WalkTest) {
  // When I am in a hurry, I walk
  Length<double> length = Length<double>::kilometer(1.0);
  // in
  Duration<double> time = Duration<double>::minute(10);

  // Then my speed is
  Velocity<double> velA = Velocity<double>::kilometerPerHour(length.kilometer()/time.hour());
  Velocity<double> velB = Velocity<double>::meterPerSecond(length.meter()/time.second());
  Velocity<double> velC = length / time;

  EXPECT_NEAR(velA.knot(), velB.knot(), 1.0e-5);
  EXPECT_NEAR(velB.meterPerSecond(), velA.meterPerSecond(), 1e-5);
}

TEST(PhysQuantTest, AngleTest) {
  Angle<double> a = Angle<double>::degree(30);
  Angle<double> b = Angle<double>::degree(60);
  Angle<double> a2 = Angle<double>::radian(3.14159265359 /6.0);
  Angle<double> b2 = Angle<double>::radian(3.14159265359 /3.0);
  EXPECT_NEAR(a.degree(), a2.degree(), 1.0e-5);
  EXPECT_NEAR(b.degree(), b2.degree(), 1.0e-5);
  EXPECT_NEAR(sin(a), 1.0/2.0, 1.0e-6);
  EXPECT_NEAR(sin(b), sqrt(3)/2.0, 1.0e-6);
  EXPECT_NEAR(sin(a)*sin(a) + cos(a)*cos(a), 1.0, 1.0e-6);
  EXPECT_NEAR(cos(a2)*cos(b2) - sin(a2)*sin(b2), cos(a.radian() + b.radian()), 1.0e-5);
}

// Try the operators
TEST(PhysQuantTest, OperatorTest) {
  Mass<double> a = Mass<double>::kilogram(30.0);
  Mass<double> b = Mass<double>::kilogram(34.0);

  // The + operator is inherited
  Mass<double> c = a + b;

  EXPECT_NEAR(c.kilogram(), 64.0, 1.0e-6);
}

TEST(PhysQuantTest, OperatorTest2) {
  Mass<double> sum = Mass<double>::lispund(0.0);
  for (int i = 0; i < 20; i++) {
    sum = sum + Mass<double>::lispund(1.0);
  }
  EXPECT_NEAR(sum.skeppund(), 1.0, 1.0e-6);
}

TEST(PhysQuantTest, OperatorTest6) {
  Mass<double> a = Mass<double>::lispund(1.0);
  Mass<double> b = a;
  EXPECT_NEAR((a - b).kilogram(), 0.0, 1.0e-6);
}

TEST(PhysQuantTest, DurationStr) {
  Duration<> d(
      Duration<>::week(3)
      + Duration<>::day(1)
      + Duration<>::hour(7)
      + Duration<>::minute(3)
      + Duration<>::second(5));
  EXPECT_EQ("3 weeks, 1 day, 7 hours, 3 minutes, 5 seconds", d.str());
}

TEST(PhysQuantTest, CastTest) {
    Mass<float> flt = Mass<float>::kilogram(3.0f);
    Mass<double> dbl = flt.cast<double>();
    EXPECT_EQ(Mass<double>::kilogram(3.0), dbl);
}

TEST(PhysQuantTest, HorizontalMotionTest) {
    HorizontalMotion<double> a(
        Velocity<double>::meterPerSecond(3),
        Velocity<double>::meterPerSecond(2));
    auto b(a);
    HorizontalMotion<double> twice(
        Velocity<double>::meterPerSecond(6),
        Velocity<double>::meterPerSecond(4));
    EXPECT_EQ(a, b);
    EXPECT_EQ(twice, HorizontalMotion<double>(a + b));
    EXPECT_EQ(a, HorizontalMotion<double>(
            a + HorizontalMotion<double>::zero()));
}

TEST(PhysQuantTest, HorizontalMotionPolarTest) {
    HorizontalMotion<double>toWestA(
        Velocity<double>::meterPerSecond(-1),
        Velocity<double>::meterPerSecond(0));
    auto toWestB = HorizontalMotion<double>::polar(
            Velocity<double>::meterPerSecond(1),
            Angle<double>::degree(270));
    EXPECT_NEAR(toWestA[0].meterPerSecond(), toWestB[0].meterPerSecond(), 1e-8);
    EXPECT_NEAR(toWestA[1].meterPerSecond(), toWestB[1].meterPerSecond(), 1e-8);

    HorizontalMotion<double>toSouthA(
        Velocity<double>::meterPerSecond(0),
        Velocity<double>::meterPerSecond(-1));
    auto toSouthB = HorizontalMotion<double>::polar(
            Velocity<double>::meterPerSecond(1),
            Angle<double>::degree(180));
    EXPECT_NEAR(toSouthA[0].meterPerSecond(), toSouthB[0].meterPerSecond(), 1e-8);
    EXPECT_NEAR(toSouthA[1].meterPerSecond(), toSouthB[1].meterPerSecond(), 1e-8);
}

TEST(PhysQuantTest, HorizontalMotionAngleNormTest) {
    for (double angle = -180 + 15; angle<180; angle += 15) {
        auto motion = HorizontalMotion<double>::polar(
                Velocity<double>::knot(10),
                Angle<double>::degree(angle));
		EXPECT_NEAR(angle, motion.angle().degree(), 1e-5)
			<< "degree: " << motion.angle().degree();
        EXPECT_NEAR(10, motion.norm().knot(), 1e-5);
    }
}

TEST(PhysQuantTest, isNaN) {
  EXPECT_FALSE(isNaN(0));
  EXPECT_FALSE(isNaN(0.0));
  EXPECT_TRUE(isNaN(std::numeric_limits<double>::signaling_NaN()));
  EXPECT_TRUE(isNaN(std::numeric_limits<float>::signaling_NaN()));
  EXPECT_TRUE(isNaN(std::numeric_limits<float>::quiet_NaN()));
}

TEST(PhysQuantTest, Literals) {
  EXPECT_NEAR((134.0_deg).degree(), 134, 1.0e-6);
}

TEST(PhysQuantTest, VectorizeNorm) {
  Vectorize<Length<double>, 2> x{3.0_m, 4.0_m};
  EXPECT_NEAR(x.norm().meter(), 5.0, 1.0e-6);
}

TEST(PhysQuantTest, AngularVelocity) {
  AngularVelocity<float> base = AngularVelocity<float>::radianPerSecond(4);
  EXPECT_NEAR(38.19718, base.rpm(), 1e-5);
  EXPECT_NEAR(4 * 180 / 3.14159265359, base.degreePerSecond(), 1e-5);
}

TEST(PhysQuantTest, TemperatureTest) {
  Temperature<double> cold = Temperature<double>::celsius(-30);
  EXPECT_NEAR(cold.kelvin(), -30 + 273.15, 1e-5);
  EXPECT_NEAR(Temperature<>::kelvin(0).celsius(), -273.15, 1e-5);
  EXPECT_NEAR(Temperature<>::celsius(0).kelvin(), 273.15, 1e-5);
  EXPECT_NEAR(Temperature<>::farenheit(123).celsius(), 50.555, 1e-3);
}

TEST(PhysQuantTest, Voltage) {
  Voltage<> u = Voltage<>::volt(12);
  Resistance<> r = Resistance<>::ohm(10);
  ElectricCurrent<> i = ElectricCurrent<>::ampere(1.2);

  EXPECT_NEAR(u.volt(), (r * i).volt(), 1e-5);
}

TEST(PhysQuantTest, Sqrt) {
  Length<double> l = sqrt(Area<>::squareMeter(10 * 10));
}
