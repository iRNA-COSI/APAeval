# OpenEBench requirements

The following schema groups need to be provided in order to register benchamrks
with [OpenEBench][oeb]:

* **Tools** and **References** (to tools): can be filled out right away from
  our [list of methods][methods].
* **Challenges**/**Metrics**: pretty much correspond to what we call
  _benchmarks_. Each challenge _can_ have multiple metrics (say you want to
  create more than one plot per benchmark). It should become clear from the
  examples and most of the required info should already be available in our
  benchmark specifications.
* **Contacts**: lists the people who have developed the challenges/metrics,
  i.e., our metrics team, as well as the contacts of method developers. The
  latter should be available in our [methods list][methods].

Please create and fill out the schemas _as best as you can_ and put them in
the appropriate subdirectory in directory [`APAeval_schemas`](APAeval_schemas).
Use the schemas in the [`examples`](examples) directory as a reference
(unfortunately the links to the schema definitions as indicated in the examples
are currently dead). To see how a complete case looks like, including the
directory structure, here is a [full set of examples][example-full].

[example-full]: <https://github.com/inab/benchmarking-data-model/tree/master/prototype-data/1.0.x/QfO>
[methods]: <https://docs.google.com/spreadsheets/d/1jvZOpF8iKzFRltrg99dgmkcHNXQVvaUla2JEjt5vyCU/edit>
[oeb]: <https://openebench.bsc.es/>
