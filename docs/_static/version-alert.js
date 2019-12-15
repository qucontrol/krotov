"use strict";

// Source:
// https://github.com/anymail/django-anymail/blob/4c443f5515d1d5269a95cb54cf75057c56a3b150/docs/_static/version-alert.js

function warnOnLatestVersion() {

  var warning = document.createElement('div');
  warning.setAttribute('class', 'admonition danger');
  warning.innerHTML = "<p class='first admonition-title'>Note</p> " +
    "<p class='last'> " +
    "The documentation for the krotov package has moved and is now available at <a href='https://qucontrol.github.io/krotov/'>https://qucontrol.github.io/krotov/</a>" +
    "</p>";

  var parent = document.querySelector('div.body')
    || document.querySelector('div.document')
    || document.body;
  parent.insertBefore(warning, parent.firstChild);
}

document.addEventListener('DOMContentLoaded', warnOnLatestVersion);
