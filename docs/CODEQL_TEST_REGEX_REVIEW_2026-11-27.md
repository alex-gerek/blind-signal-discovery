# CodeQL test-regex review and disposition (2026-11-27)

Owner approved the previously reported findings and closure before any disposition change. Alerts1 and2 remained open until this approval.

- Rule:js/bad-tag-filter, security severityHIGH, CodeQL classificationtest.
- Reviewed main:5e19ae4e6fb9ea177738fa4bae8310c11a2633e6.
- Alert1:tests/run_browser_wav_resampler.js:5. Fixed repositoryHTML is read/extracted/executed before stdinJSON is read for numeric resampling input.
- Alert2:tests/test_browser_input_validation.js:9. Fixed repositoryHTML is extracted for a mocked-DOM Node harness; uploaded/testCSV/WAV does not flow into the extraction regex.
- Source HTML and test files were fetched from the exact main revision and matched local source. Pages servesmain:/docs, with no loading of these Node test harnesses; live HTML matched the reviewed source.
- Reproduced:uppercaseSCRIPT is not matched and extraction throwsTypeError. The regexp is brittle against changes to the trusted fixture format, but it is not an application HTML sanitizer or an untrusted-input execution boundary.
- Earlier read-only investigation passed16Node scenarios and3SciPy/pytest reference cases; source unchanged. These checks were not rerun merely to change alert metadata.
- Conclusion:security-context false positives with acknowledged test-format fragility; no claim that the entire application is free of XSS.

On2026-11-27 both alerts were changed to dismissed with reasonfalse positive and this comment:

Test-only extraction from fixed repository HTML, not sanitization. No uploaded/stdin HTML reaches this regex; Pages does not load the Node harness. Uppercase SCRIPT fails extraction. Main 5e19ae4 reviewed; 16 Node and 3 SciPy checks passed.

[Alert1](https://github.com/alex-gerek/blind-signal-discovery/security/code-scanning/1), [alert2](https://github.com/alex-gerek/blind-signal-discovery/security/code-scanning/2). A first PATCH attempt using a longer comment failed before disposition was changed; the reason text was not retained. A subsequent240-character comment succeeded for both alerts. No code, detector, test or GitHub security setting was changed. Optional fixture-extractor cleanup remains a maintenance task, not an urgent application-security repair.
