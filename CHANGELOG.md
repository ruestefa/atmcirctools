# AtmCircTools: Changelog

## v0.?.? (2022-??-??)

- Update `atmcirclib` to `v0.4.0` and update envs
- Hook project up to copier template [`python-project`](https://git.iac.ethz.ch/atmcirc/templates/python-project) and update
- Remove leftovers from former scikit-build project structure
- Implement new command `act check-ncfiles` to sanity-check NetCDF files (by opening them)
- Remove dependency `h5netcdf`; not used directly, but caused problems with editable pip-install in conda env (pip identified version as 0.0.0, not 1.0.0)

## v0.2.2 (2022-06-08)

- Add missing requirement `toml` and indirect requirement `scikit-build` (for `create-recipe`) and update envs
- Remove archived conda build recipes; only retain latest as `recipe/meta.yaml`
- Change default output file of `create-recipe` from `recipe/v{version}/meta.yaml` to `recipe/meta.yaml`

## ~v0.2.1 (2022-06-08)~

- Incomplete; superseded by v0.2.2

## v0.2.0 (2022-05-31)

- Add script `atmcirclib.bin.call_graph` as command `act call-graph`
- Add script `atmcirclib.bin.create_recipe` as command `act create-recipe`
- Remove command line tool `act-intp-lvl`

## v0.1.0 (2022-04-25)

- Write command line tool `act-intp-lvl` to interpolate a 3D field to vertical levels (e.g., PV to TH surfaces)
- Write umbrella executable `act` (short for `AtmCircTools`) that provides `act-intp-lvl` as `act intp-lvl`
