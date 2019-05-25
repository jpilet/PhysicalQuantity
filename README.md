# Introduction

Physical Quantity is a single-header c++ library providing compile-time SI unit consistency.

Physical Quantity implements most of the International System of Units (SI),
with a few exceptions (angles, for example).  Once you start using it, your
compiler will take care of making sure all units are converted properly, at
compile time.

Here's an example:

```
  // When I am in a hurry, I walk
  Length<double> length = Length<double>::kilometer(1.0);
  // in
  Duration<double> time = Duration<double>::minute(10);

  // Then my speed is
  Velocity<double> velC = length / time;

  std::cout << "My walking speed: " << velC.meterPerSecond() << "m/s, "
   << velC.kilometerPerHour() << "km/h, "
   << velC.knots() << " knots.\n";
```

Basically, the library provides template types to wrap the containing values
and do the proper conversions.

When using the variable, there is no need to worry about its unit: it is enough
to know its nature (a length, a duration, a voltage,...).

# Other libraries solving a similar problem

- boost units: more complete and flexible, more complicated to use.
- [unit lite](https://github.com/pierreblavy2/unit_lite): allows you to define your own units, but less convenient for SI.
- [Poco Units](https://pocoproject.org/docs/Poco.Util.Units.html).


# What quantities/units are supported?

Most of SI and part of the imperial system are supported, with both British (meter) and American (meter) spellings.

See the [full list](./UNITS.md)

# What operator is supported?
- multiplication with another quantity ```*```
- division ```/```
- addition ```+```
- subtraction  ```-```
- multiplication with a scalar (with ```quantity.scaled(scalar)``` or ```scalar * quantity```)
- sqrt
- fabs
- Trigonometry: cos, sin, tan
- isNaN, isFinite

# How to use Physical Quantity in my project?

You need just a single file: PhysicalQuantity.h. Make sure you can include it,
either by copying it in your project or add the containing folder to the
include directories.

