# APAeval

Welcome to the [APAeval][apa-eval] GitHub repository.

APAeval is a community effort to evaluate computational methods for the
detection and quantification of poly(A) sites and the estimation of their
differential usage across RNA-seq samples.

[![logo][apa-eval-logo]][apa-eval]

## What is APAeval?

APAeval consists of three challenges, each consisting of a set of benchmarks
for tools that:

1. **Identify** polyadenylation sites
2. **Quantify** polyadenylation sites
3. Calculate **differential usage** of polyadenylation sites

For more info, please refer to our [landing page][apa-eval].

## How to get involved?

If you would like to contribute to APAeval, the first things you would need to
do are:

- [Register][apa-eval] for one or more of the APAeval teams
- Please use [this form][form-service-accounts] to provide us with your user
  names/handles of the various service accounts we are using (see
  [below](#how-do-we-work)) so that we can add you to the corresponding
  repositories/organizations
- Wait for us to reach out to you with an invitation to our Slack space, then
  use the invitation to sign up and post a short intro message in the
  [`#general`][slack-general] channel

## How do we work?

To account for everyone's different agendas and time zones, the hackathon is
organized such that contributors can work, as much as possible, in their own
time.

### Communication

#### Chat

We are making use of [**Slack**][slack] (see [above](#how-to-get-involved) to
see how you can become a member) for asynchronous communication. Please use the
most appropriate channels for discussions/questions etc.:

- [`#general`][slack-general]: Introduce yourself and find pointers to get you
  started. All APAeval-wide announcements will be put here!
- [`#admin`][slack-admin]: Ask questions about the general organization of
  APAeval.
- [`#tech-support`][slack-tech-support]: Ask questions about the
  technical infrastructure and relevant software, e.g., AWS, GitHub, Nextflow,
  Snakemake.
- [`#challenge-differential-usage`][slack-diff-usage]: Discuss the poly(A) site
  differential usage challenge.
- [`#challenge-identification`][slack-identification]: Discuss the poly(A) site
  identification challenge.
- [`#challenge-quantification`][slack-quantification]: Discuss the poly(A) site
  quantification challenge.
- `#tool-{METHOD_NAME}`: Discussions channels for each indidvidual method.
- [`#random`][slack-random]: Post anything that doesn't fit into any of the
  other channels.

#### Video calls

Despite the event taking place mostly asynchrounously, we do have a few video
calls to increase the feeling of collaboration. In particular, we hae a kickoff
and a wrap-up meeting, three reporting calls and the presentation of
APAeval during the RNA Meeting. Please join meetings if you can, we have tried
to shuffle the times a bit, so that most people should be able to attend at
least some of the calls, irrespective of your location.

This calendar contains all video call events, including the necessary login
info, and we would like to kindly ask you to subscribed to it:

- Calendar ID: `59bboug9agv30v32r6bvaofdo4@group.calendar.google.com`
- Public address:
  `https://calendar.google.com/calendar/ical/59bboug9agv30v32r6bvaofdo4%40group.calendar.google.com/public/basic.ics`

> Please do not download the ICS file and then import it, as any updates to the
> calendar will not be synced. Refer to your calendar application's help pages
> if you do not know how to _subscribe_ to a calendar.

With the exception of presenting APAeval during the RNA Meeting, all video
calls will take place in the following [**Zoom**][zoom] room:

- [Direct link][vc-direct-link]
- Meeting ID: `656 9429 1427`
- Passcode: `APAeval`

There is also a [meeting agenda][agenda-vc].

> For more lively meetings, participants are encouraged to switch on their
> cameras. But please mute your microphones if you are not currently speaking.

### Social coding

We are making extensive use of [**GitHub**][gh]'s project management resources
to allow contributors to work independently on individual, largely
self-contained, issues. There are three kanban project boards, one for each
challenge:

- [Identification challenge][gh-project-identification]
- [Quantification challenge][gh-project-quantification]
- [Differential usage challenge][gh-project-diff-usage]

The idea is that people assign themselves to open issues (i.e., those issues
that are not yet assigned to someone else). Note that in order to be able to
do so, you will need to be a member of this GitHub repository (see
[above](#how-to-get-involved) to see how you can become a member). Once you
have assigned yourself, you can move/drag the issue from the "To do" to the "In
progress" column of the kanban board. It is not necessary to 

When working on an issue, please start by cloning (preferred) or forking) the
repository. Then create a feature branch and implement your code/changes. Once
you are happy with them, please create a pull request against the `main`
branch, making sure to fill in the provided template (in particular, please
refer to the original issue you are addressing with this pull request) and to
assign two reviewers. This workflow ensures collaborative coding and is
sometimes referred to as [GitHub flow][gh-flow]. If you are not familiar with
Git, GitHub or the GitHub flow, there are many useful tutorials online.

### Storage & compute infrastructure

TBA

[apa-eval]: <https://irnacosi.org/2021/01/04/rna-society-2021-apaeval-challenge/> 
[apa-eval-logo]: images/logo.png
[apa-eval-members]: <https://docs.google.com/document/d/1G7u-WQ6C-I_sXZ-15CIBw2iNgw6jkTNo7hnRTjci_b4/edit#heading=h.tarrapa8v8n6>
[contact]: <mailto:apaeval@irnacosi.org>
[docker-hub]: <https://hub.docker.com>
[docker-hub-org]: <https://hub.docker.com/orgs/apaeval>
[form-service-accounts]: <https://forms.gle/eKCHe5GWtvGrriek8>
[gh]: <http://github.com/>
[gh-flow]: <https://guides.github.com/introduction/flow/>
[gh-join]: <https://github.com/join>
[gh-project-identification]: <https://github.com/iRNA-COSI/APAeval/projects/1>
[gh-project-quantification]: <https://github.com/iRNA-COSI/APAeval/projects/2>
[gh-project-diff-usage]: <https://github.com/iRNA-COSI/APAeval/projects/4>
[slack]: <http://slack.com/>
[slack-admin]: <https://apaeval.slack.com/archives/C01PEJQEUMT>
[slack-diff-usage]: <https://apaeval.slack.com/archives/C022DBAFSE7>
[slack-identification]: <https://apaeval.slack.com/archives/C0232UJQREU>
[slack-quantification]: <https://apaeval.slack.com/archives/C022A2PM4GM>
[slack-general]: <https://apaeval.slack.com/archives/C01PHLQKNH0>
[slack-random]: <https://apaeval.slack.com/archives/C01Q7FMRJ3A>
[slack-tech-support]: <https://apaeval.slack.com/archives/C022RNSAUV7>
[agenda-vc]: <https://docs.google.com/document/d/1Cl3xq7_uwApAYxUbzeVSBsRfGUmtRc0jSnZ3yrWM3ks/edit#>
[vc-direct-link]: <https://unibas.zoom.us/j/65694291427?pwd=QUMyMjQ2SSt2eS9iZW50YVZCOC8wQT09>
[zoom]: <https://zoom.us/>
