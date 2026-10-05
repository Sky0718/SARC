from fractions import Fraction

def decode_rationals(value):
    if (isinstance(value, list)):
        return [decode_rationals(item) for (item) in (value)]
    if (isinstance(value, dict)):
        fields = set(value)
        if (fields in (
            {"numerator", "denominator"},
            {"numerator", "denominator", "value"},
        )):
            numerator, denominator = value["numerator"], value["denominator"]
            if (
                type(numerator) is not int
                or type(denominator) is not int
                or denominator <= 0
            ):
                raise ValueError("Invalid exact rational carrier")
            rational = Fraction(numerator, denominator)
            if ("value" in value and (
                not isinstance(value["value"], (int, float))
                or isinstance(value["value"], bool)
                or float(rational) != value["value"]
            )):
                raise ValueError(
                    "Display value differs from its exact rational carrier"
                )
            return rational
        return {key: decode_rationals(item) for ((key, item)) in (value.items())}
    return value
