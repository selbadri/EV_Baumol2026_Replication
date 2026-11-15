# Jordà-Schularick-Taylor Macrohistory Database

Used only by `code/07c_markups.do`, to build the required real return on wealth
for the machine-producer profit shares (Appendix F).

## Source file

| File | Source | Content |
|---|---|---|
| `JSTdatasetR6.dta` | Jordà-Schularick-Taylor Macrohistory Database, release 6 (www.macrohistory.net), downloaded 1 October 2026 | Annual macro-financial series for 18 economies, 1870-2020 |

Variables used (United States only): `eq_tr`, `housing_tr`, `bond_tr`, `bill_rate`,
`capital_tr` (returns), `eq_dp` (dividend-price ratio), `housing_rent_yd` (rental yield),
`ltrate`, `cpi`, `rgdpbarro`, `pop`, `hpnom`.

## Use

The required return is a yield-based expected real return on wealth: dividend yield plus
mean real GDP growth, rental yield plus mean real house-price growth (both means over
1950-2020), and the long-term and bill rates less expected inflation (mean CPI inflation of
the previous five years). The four components are weighted by the coefficients of a
regression of JST's realized return on wealth (`capital_tr`) on the four asset returns over
1950-2020, without a constant, rescaled to sum to one (0.248, 0.486, 0.084, 0.181). The
series ends in 2020 and is held at its 2020 value for 2021-2023. Mean 1950-2023: 4.7
percent.

## License and citation

The database is distributed under a Creative Commons BY-NC-SA 4.0 license. Cite:

- Jordà, Òscar, Moritz Schularick, and Alan M. Taylor. 2017. "Macrofinancial History and
  the New Business Cycle Facts." NBER Macroeconomics Annual 31: 213-263.
- Jordà, Òscar, Katharina Knoll, Dmitry Kuvshinov, Moritz Schularick, and Alan M. Taylor.
  2019. "The Rate of Return on Everything, 1870-2015." Quarterly Journal of Economics 134(3):
  1225-1298.
