---
title: Thing 1 topic
permalink: thing1.html
sidebar: generic
product: Generic
---

Create a file called samplelist.yml in _data.

Put something like this in it:

    toc:
      - title: Group 1
        subfolderitems:
          - page: Thing 1
            url: /thing1.html
          - page: Thing 2
            url: /thing2.html
          - page: Thing 3
            url: /thing3.html
      - title: Group 2
        subfolderitems:
          - page: Piece 1
            url: /piece1.html
          - page: Piece 2
            url: /piece2.html
          - page: Piece 3
            url: /piece3.html
      - title: Group 3
        subfolderitems:
          - page: Widget 1
            url: /widget1.html
          - page: Widget 2
            url: /widget2.html
          - page: Widget 3
            url: /widget3.html

Now loop through it like this:

    <div class="result">
    {% for item in site.data.samplelist.toc %}
    <h3>{{item.title}}</h3>
    <ul>
    {% for entry in item.subfolderitems %}
    <li><a href="{{entry.url}}">{{entry.page}}</a></li>
    {% endfor %}
    </ul>
    {% endfor %}
    </div>
