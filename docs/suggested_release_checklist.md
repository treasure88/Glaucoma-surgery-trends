# Suggested public release checklist

Before creating the GitHub release and the Zenodo DOI, confirm the following:

1. Remove all local machine paths (for example, Desktop paths or drive letters).
2. Remove all internal database connection details, server names, usernames, or passwords.
3. Remove any intermediate exports that contain patient-level rows.
4. Confirm that the codebook contains variable definitions only and no patient-level values.
5. Confirm that all plots, example files, and screenshots are free of identifiable information.
6. Replace GitHub username and DOI placeholders in:
   - README.md
   - CITATION.cff
   - zenodo.json
   - manuscript Data/Code Availability statement
7. Create a tagged GitHub release (for example, v1.0.0).
8. After Zenodo mints the DOI, update the repository badge and manuscript text.
