# Change Log

## [1.0.6] 2021-05-16
### Dependencies Update

- Freeze used versions in `requirements.txt`
    - jinja2 = 2.11.3

## [1.0.5] 2021-03-18
### Improvements

- Freeze used versions in `requirements.txt`
    - flask_sqlalchemy = 2.4.4
    - sqlalchemy = 1.3.23

## [1.0.4] 2021-01-20
### Improvements

- Remove `shutdown` route from Base Blueprint - Updated file(s):
    - `app\base\routes.py`

## [1.0.3] 2021-01-01
### Improvements 

- 2021-01-01 - Registration
    - Hide form on success

## [1.0.2] 2020-12-29
### Improvements & Bug Fixes

- 2020-12-29 - `.env` usage, add logging
    - Updated files(s): run.py

- 2020-08-20 - Added get_segment() helper that detects the current page
    - Updated files(s): app/home/routes.py

- 2020-06-22 - Guard Flask links with quotes
    - Sample href="{{ url_for('base_blueprint.login') }}"
    - Impacted files: login.html, register.html, sidebar.html

- 2020-06-22 - Added HEROKU support. Impacted files:
    - runtime.txt - Bump the Python version to 3.6.10
    - README added new section for HEROKU deployment

## [1.0.1] 2020-05-30
### Improvements & Bug Fixes

- Patch #Bug - Return a 403 Error for unauthorized access
- Update Licensing information
- Add CHANGELOG.md to track all changes

## [1.0.0] 2020-02-07
### Initial Release
