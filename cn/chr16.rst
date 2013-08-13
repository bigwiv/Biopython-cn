


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>Biopython-cn/chr16.rst at master · julyfire/Biopython-cn</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="https://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <meta property="og:image" content="https://github.global.ssl.fastly.net/images/modules/logos_page/Octocat.png">
    <meta name="hostname" content="fe3.rs.github.com">
    <meta name="ruby" content="ruby 1.9.3p194-tcs-github-tcmalloc (2012-05-25, TCS patched 2012-05-27, GitHub v1.0.32) [x86_64-linux]">
    <link rel="assets" href="https://github.global.ssl.fastly.net/">
    <link rel="xhr-socket" href="/_sockets" />
    
    


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" /><meta content="5134276" name="octolytics-actor-id" /><meta content="julyfire" name="octolytics-actor-login" /><meta content="18bef267938255a85af25df4cca3658b03e9ead17db2ef6f4e58ddea04eda9f1" name="octolytics-actor-hash" />

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="MaNXeYx+GszEHT/wsFpw/ebGbHMr5pbNULyIthkFlSY=" name="csrf-token" />

    <link href="https://github.global.ssl.fastly.net/assets/github-697b506a0b29d9891d9591c88950d8c4ab0c7c0b.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://github.global.ssl.fastly.net/assets/github2-d3cdff01428383b2670dce9b781434be7c9d37e8.css" media="all" rel="stylesheet" type="text/css" />
    


      <script src="https://github.global.ssl.fastly.net/assets/frameworks-eae23340ab2a6ba722166712e699c87aaf60ad8f.js" type="text/javascript"></script>
      <script src="https://github.global.ssl.fastly.net/assets/github-a4de942c73c48e278373120f94de228c91dde62a.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="a7b574fbd4948374ae3230fda6444452">

        <link data-pjax-transient rel='permalink' href='/julyfire/Biopython-cn/blob/c33a446bfbd1d8d85f49f8c5455a48f1c644c718/chr16.rst'>
  <meta property="og:title" content="Biopython-cn"/>
  <meta property="og:type" content="githubog:gitrepository"/>
  <meta property="og:url" content="https://github.com/julyfire/Biopython-cn"/>
  <meta property="og:image" content="https://github.global.ssl.fastly.net/images/gravatars/gravatar-user-420.png"/>
  <meta property="og:site_name" content="GitHub"/>
  <meta property="og:description" content="Biopython-cn - Biopython 中文翻译"/>

  <meta name="description" content="Biopython-cn - Biopython 中文翻译" />

  <meta content="5134276" name="octolytics-dimension-user_id" /><meta content="julyfire" name="octolytics-dimension-user_login" /><meta content="12048834" name="octolytics-dimension-repository_id" /><meta content="julyfire/Biopython-cn" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="true" name="octolytics-dimension-repository_is_fork" /><meta content="10706983" name="octolytics-dimension-repository_parent_id" /><meta content="luwening/Biopython-cn" name="octolytics-dimension-repository_parent_nwo" /><meta content="10706983" name="octolytics-dimension-repository_network_root_id" /><meta content="luwening/Biopython-cn" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/julyfire/Biopython-cn/commits/master.atom" rel="alternate" title="Recent Commits to Biopython-cn:master" type="application/atom+xml" />

  </head>


  <body class="logged_in page-blob windows vis-public fork env-production ">

    <div class="wrapper">
      
      
      


      <div class="header header-logged-in true">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/">
  <span class="mega-octicon octicon-mark-github"></span>
</a>

    <div class="divider-vertical"></div>

      <a href="/notifications" class="notification-indicator tooltipped downwards" title="You have no unread notifications">
    <span class="mail-status all-read"></span>
  </a>
  <div class="divider-vertical"></div>


      <div class="command-bar js-command-bar  in-repository">
          <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<input type="text" data-hotkey="/ s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    data-username="julyfire"
      data-repo="julyfire/Biopython-cn"
      data-branch="master"
      data-sha="26b073e57c8eb98a96b744e3d2a2aa84170d5363"
  >

    <input type="hidden" name="nwo" value="julyfire/Biopython-cn" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
        <ul class="top-nav">
            <li class="explore"><a href="/explore">Explore</a></li>
            <li><a href="https://gist.github.com">Gist</a></li>
            <li><a href="/blog">Blog</a></li>
          <li><a href="https://help.github.com">Help</a></li>
        </ul>
      </div>

    

  

    <ul id="user-links">
      <li>
        <a href="/julyfire" class="name">
          <img height="20" src="https://secure.gravatar.com/avatar/588d555de56976d4893c5dd49a4e260a?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png" width="20" /> julyfire
        </a>
      </li>

        <li>
          <a href="/new" id="new_repo" class="tooltipped downwards" title="Create a new repo" aria-label="Create a new repo">
            <span class="octicon octicon-repo-create"></span>
          </a>
        </li>

        <li>
          <a href="/settings/profile" id="account_settings"
            class="tooltipped downwards"
            aria-label="Account settings "
            title="Account settings ">
            <span class="octicon octicon-tools"></span>
          </a>
        </li>
        <li>
          <a class="tooltipped downwards" href="/logout" data-method="post" id="logout" title="Sign out" aria-label="Sign out">
            <span class="octicon octicon-log-out"></span>
          </a>
        </li>

    </ul>


<div class="js-new-dropdown-contents hidden">
  

<ul class="dropdown-menu">
  <li>
    <a href="/new"><span class="octicon octicon-repo-create"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new"><span class="octicon octicon-organization"></span> New organization</a>
  </li>



    <li class="section-title">
      <span title="julyfire/Biopython-cn">This repository</span>
    </li>
    <li>
      <a href="/julyfire/Biopython-cn/issues/new"><span class="octicon octicon-issue-opened"></span> New issue</a>
    </li>
      <li>
        <a href="/julyfire/Biopython-cn/settings/collaboration"><span class="octicon octicon-person-add"></span> New collaborator</a>
      </li>
</ul>

</div>


    
  </div>
</div>

      

      




          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">

    <li class="subscription">
      <form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="MaNXeYx+GszEHT/wsFpw/ebGbHMr5pbNULyIthkFlSY=" /></div>  <input id="repository_id" name="repository_id" type="hidden" value="12048834" />

    <div class="select-menu js-menu-container js-select-menu">
        <a class="social-count js-social-count" href="/julyfire/Biopython-cn/watchers">
          1
        </a>
      <span class="minibutton select-menu-button with-count js-menu-target">
        <span class="js-select-button">
          <span class="octicon octicon-eye-unwatch"></span>
          Unwatch
        </span>
      </span>

      <div class="select-menu-modal-holder">
        <div class="select-menu-modal subscription-menu-modal js-menu-content">
          <div class="select-menu-header">
            <span class="select-menu-title">Notification status</span>
            <span class="octicon octicon-remove-close js-menu-close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-list js-navigation-container">

            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_included" name="do" type="radio" value="included" />
                <h4>Not watching</h4>
                <span class="description">You only receive notifications for discussions in which you participate or are @mentioned.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-watch"></span>
                  Watch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input checked="checked" id="do_subscribed" name="do" type="radio" value="subscribed" />
                <h4>Watching</h4>
                <span class="description">You receive notifications for all discussions in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-unwatch"></span>
                  Unwatch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_ignore" name="do" type="radio" value="ignore" />
                <h4>Ignoring</h4>
                <span class="description">You do not receive any notifications for discussions in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-mute"></span>
                  Stop ignoring
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

</form>
    </li>

  <li>
  
<div class="js-toggler-container js-social-container starring-container ">
  <a href="/julyfire/Biopython-cn/unstar" class="minibutton with-count js-toggler-target star-button starred upwards" title="Unstar this repo" data-remote="true" data-method="post" rel="nofollow">
    <span class="octicon octicon-star-delete"></span><span class="text">Unstar</span>
  </a>
  <a href="/julyfire/Biopython-cn/star" class="minibutton with-count js-toggler-target star-button unstarred upwards " title="Star this repo" data-remote="true" data-method="post" rel="nofollow">
    <span class="octicon octicon-star"></span><span class="text">Star</span>
  </a>
  <a class="social-count js-social-count" href="/julyfire/Biopython-cn/stargazers">0</a>
</div>

  </li>


        <li>
          <a href="/julyfire/Biopython-cn/fork" class="minibutton with-count js-toggler-target fork-button lighter upwards" title="Fork this repo" rel="facebox nofollow">
            <span class="octicon octicon-git-branch-create"></span><span class="text">Fork</span>
          </a>
          <a href="/julyfire/Biopython-cn/network" class="social-count">9</a>
        </li>


</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo-forked"></span>
          <span class="author">
            <a href="/julyfire" class="url fn" itemprop="url" rel="author"><span itemprop="title">julyfire</span></a></span
          ><span class="repohead-name-divider">/</span><strong
          ><a href="/julyfire/Biopython-cn" class="js-current-repository js-repo-home-link">Biopython-cn</a></strong>

          <span class="page-context-loader">
            <img alt="Octocat-spinner-32" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

            <span class="fork-flag">
              <span class="text">forked from <a href="/luwening/Biopython-cn">luwening/Biopython-cn</a></span>
            </span>
        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">

      <div class="repository-with-sidebar repo-container ">

        <div class="repository-sidebar">
            

<div class="repo-nav repo-nav-full js-repository-container-pjax js-octicon-loaders">
  <div class="repo-nav-contents">
    <ul class="repo-menu">
      <li class="tooltipped leftwards" title="Code">
        <a href="/julyfire/Biopython-cn" aria-label="Code" class="js-selected-navigation-item selected" data-gotokey="c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /julyfire/Biopython-cn">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


      <li class="tooltipped leftwards" title="Pull Requests"><a href="/julyfire/Biopython-cn/pulls" aria-label="Pull Requests" class="js-selected-navigation-item js-disable-pjax" data-gotokey="p" data-selected-links="repo_pulls /julyfire/Biopython-cn/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped leftwards" title="Wiki">
          <a href="/julyfire/Biopython-cn/wiki" aria-label="Wiki" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_wiki /julyfire/Biopython-cn/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="repo-menu-separator"></div>
    <ul class="repo-menu">

      <li class="tooltipped leftwards" title="Pulse">
        <a href="/julyfire/Biopython-cn/pulse" aria-label="Pulse" class="js-selected-navigation-item " data-pjax="true" data-selected-links="pulse /julyfire/Biopython-cn/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Graphs">
        <a href="/julyfire/Biopython-cn/graphs" aria-label="Graphs" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_graphs repo_contributors /julyfire/Biopython-cn/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Network">
        <a href="/julyfire/Biopython-cn/network" aria-label="Network" class="js-selected-navigation-item js-disable-pjax" data-selected-links="repo_network /julyfire/Biopython-cn/network">
          <span class="octicon octicon-git-branch"></span> <span class="full-word">Network</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

    </ul>

      <div class="repo-menu-separator"></div>
      <ul class="repo-menu">
        <li class="tooltipped leftwards" title="Settings">
          <a href="/julyfire/Biopython-cn/settings" data-pjax aria-label="Settings">
            <span class="octicon octicon-tools"></span> <span class="full-word">Settings</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </a>
        </li>
      </ul>
  </div>
</div>

            <div class="only-with-full-nav">
              

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=push">
  <h3><strong>HTTPS</strong> clone URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/julyfire/Biopython-cn.git" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/julyfire/Biopython-cn.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>

  

<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=push">
  <h3><strong>SSH</strong> clone URL</h3>

  <input type="text" class="clone js-url-field"
         value="git@github.com:julyfire/Biopython-cn.git" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="git@github.com:julyfire/Biopython-cn.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=push">
  <h3><strong>Subversion</strong> checkout URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/julyfire/Biopython-cn" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/julyfire/Biopython-cn" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>



<p class="clone-options">You can clone with
    <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
    <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>,
    <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>,
  and <a href="https://help.github.com/articles/which-remote-url-should-i-use">other methods.</a>
</p>


  <a href="http://windows.github.com" class="minibutton sidebar-button">
    <span class="octicon octicon-device-desktop"></span>
    Clone in Desktop
  </a>

                <a href="/julyfire/Biopython-cn/archive/master.zip"
                   class="minibutton sidebar-button"
                   title="Download this repository as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
            </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<!-- blob contrib key: blob_contributors:v21:145de6b986b342e371559236108af3c8 -->
<!-- blob contrib frag key: views10/v8/blob_contributors:v21:145de6b986b342e371559236108af3c8 -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<a href="/julyfire/Biopython-cn/find/master" data-pjax data-hotkey="t" style="display:none">Show File Finder</a>

<div class="file-navigation">
  


<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-remove-close js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Find or create a branch…">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/julyfire/Biopython-cn/blob/master/chr16.rst" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="master" data-skip-pjax="true" rel="nofollow" title="master">master</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <form accept-charset="UTF-8" action="/julyfire/Biopython-cn/branches" class="js-create-branch select-menu-item select-menu-new-item-form js-navigation-item js-new-item-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="MaNXeYx+GszEHT/wsFpw/ebGbHMr5pbNULyIthkFlSY=" /></div>
            <span class="octicon octicon-git-branch-create select-menu-item-icon"></span>
            <div class="select-menu-item-text">
              <h4>Create branch: <span class="js-new-item-name"></span></h4>
              <span class="description">from ‘master’</span>
            </div>
            <input type="hidden" name="name" id="name" class="js-new-item-value">
            <input type="hidden" name="branch" id="branch" value="master" />
            <input type="hidden" name="path" id="branch" value="chr16.rst" />
          </form> <!-- /.select-menu-item -->

      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/julyfire/Biopython-cn" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">Biopython-cn</span></a></span></span><span class="separator"> / </span><strong class="final-path">chr16.rst</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="chr16.rst" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


  <div class="commit commit-loader file-history-tease js-deferred-content" data-url="/julyfire/Biopython-cn/contributors/master/chr16.rst">
    Fetching contributors…

    <div class="participation">
      <p class="loader-loading"><img alt="Octocat-spinner-32-eaf2f5" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32-EAF2F5.gif" width="16" /></p>
      <p class="loader-error">Cannot retrieve contributors at this time</p>
    </div>
  </div>

<div id="files" class="bubble">
  <div class="file">
    <div class="meta">
      <div class="info">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
          <span>638 lines (517 sloc)</span>
        <span>28.503 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
                <a class="minibutton"
                   href="/julyfire/Biopython-cn/edit/master/chr16.rst"
                   data-method="post" rel="nofollow" data-hotkey="e">Edit</a>
          <a href="/julyfire/Biopython-cn/raw/master/chr16.rst" class="button minibutton " id="raw-url">Raw</a>
            <a href="/julyfire/Biopython-cn/blame/master/chr16.rst" class="button minibutton ">Blame</a>
          <a href="/julyfire/Biopython-cn/commits/master/chr16.rst" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
            <a class="minibutton danger empty-icon tooltipped downwards"
               href="/julyfire/Biopython-cn/delete/master/chr16.rst"
               title="" data-method="post" rel="nofollow">
            Delete
          </a>
      </div><!-- /.actions -->

    </div>
      
  <div id="readme" class="blob instapaper_body">
    <article class="markdown-body entry-content" itemprop="mainContentOfPage">第十六章 监督学习方法
=======================================

注意本章介绍的所有监督学习方法都需要先安装Numerical Python （numpy）。

16.1 Logistic 回归模型
-----------------------------------

16.1.1 背景和目的
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Logistic回归是一种监督学习方法，通过若干预测变量*x*\ :sub:`*i*`的加权和尝试将样本划分为*K*个不同类别。Logistic回归模型可用来计算预测变量的权重βi。在Biopython中，logistic回归模型目前只实现了二类别（*K*=2）分类，而预测变量的数量没有限制。

作为一个例子，我们试着预测细菌中的操纵子结构。一个操纵子是在一条DNA链上许多相邻基因组成的一个集合，可以被共同转录为一条mRNA分子，这条mRNA分子经翻译后产生多个不同的蛋白质。我们将以枯草芽孢杆菌的操纵子数据进行说明，它的一个操纵子平均包含2.4个基因。

首先为了理解细菌的基因调节，我们需要知道其操纵子的结构。枯草芽孢杆菌大约10%的基因操纵子结构已经通过实验获知。剩下的90%的基因操纵子结构可以通过一种监督学习方法来预测。

在这种方法中，我们需要选择某些与操纵子结构有关的容易计算的预测变量*x*\ :sub:`*i*`。例如可以选择基因间碱基对距离来来作为其中一个预测变量。同一个操纵子中的相邻基因往往距离相对较近，而位于不同操纵子的相邻基因间通常具有更大的空间来容纳启动子和终止子序列。另一个预测变量可以基于基因表达量度。根据操纵子的定义，属于同一个操纵子的基因有相同的基因表达轮廓，而不同操纵子的两个基因的表达轮廓也不相同。在实际操作中，由于存在测量误差，对相同操纵子的基因表达轮廓的测量不会完全一致。为了测量基因表达轮廓的相似性，我们假设测量误差服从正态分布，然后计算对应的对数似然分。

现在我们有了两个预测变量，可以据此预测在同一条DNA链上两个相邻基因是否属于相同的操纵子：
-  *x*\ :sub:`1`：两基因间的碱基对数
-  *x*\ :sub:`2`：两基因表达轮廓的相似度

在logistic回归模型中，我们使用这两个预测变量的加权和来计算一个联合得分*S*：

+-----------------------------------------------------------------------------------+
| *S* = β:sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`.     (16.1)   |
+-----------------------------------------------------------------------------------+

根据下面两组示例基因，logistic回归模型对参数β :sub:`0`, β :sub:`1`, β :sub:`2`给出合适的值：
-  OP: 相邻基因，相同DNA链，属于相同操纵子
-  NOP: 相邻基因，相同DNA链，属于不同操纵子

在logistic回归模型中，属于某个类别的概率依赖于通过logistic函数得出的分数。对于这两类OP和NOP，相应概率可如下表述

     

Pr(\ *OP*\ \|\ *x*\ :sub:`1`, \ *x*\ :sub:`2`)

 =

 

+--------------------------------------------------------------------------+
| exp(β\ :sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`)     |
+--------------------------------------------------------------------------+
+--------------------------------------------------------------------------+
| 1+exp(β\ :sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`)   |
+--------------------------------------------------------------------------+

   

    (16.2)

Pr(\ *NOP*\ \|\ *x*\ :sub:`1`, \ *x*\ :sub:`2`)

 =

 

+--------------------------------------------------------------------------+
| 1                                                                        |
+--------------------------------------------------------------------------+
+--------------------------------------------------------------------------+
| 1+exp(β\ :sub:`0` + β:sub:`1` *x*\ :sub:`1` + β:sub:`2` *x*\ :sub:`2`)   |
+--------------------------------------------------------------------------+

   

    (16.3)

使用一组已知是否属于相同操纵子（OP类别）或不同操纵子（NOP类别）的基因对，通过最大化相应概率函数(`16.2 `__) and (`16.3 `__)的对数似然值，我们可以计算权重 β\ :sub:`0`, β\ :sub:`1`, β\ :sub:`2`.

16.1.2 训练logistic回归模型
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    --------------

    +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | 表16.1： 已知类别(OP or NOP)的相邻基因对.如果两个基因相重叠，其基因间距离为负值.   |
    +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | Gene pair           | Intergene distance (*x*\ :sub:`1`)   | Gene expression score (*x*\ :sub:`2`)   | Class   |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *cotJA* — *cotJB*   | -53                                  | -200.78                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yesK* — *yesL*     | 117                                  | -267.14                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplA* — *lplB*     | 57                                   | -163.47                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplB* — *lplC*     | 16                                   | -190.30                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplC* — *lplD*     | 11                                   | -220.94                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lplD* — *yetF*     | 85                                   | -193.94                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfmT* — *yfmS*     | 16                                   | -182.71                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfmF* — *yfmE*     | 15                                   | -180.41                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *citS* — *citT*     | -26                                  | -181.73                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *citM* — *yflN*     | 58                                   | -259.87                                 | OP      |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfiI* — *yfiJ*     | 126                                  | -414.53                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *lipB* — *yfiQ*     | 191                                  | -249.57                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfiU* — *yfiV*     | 113                                  | -265.28                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yfhH* — *yfhI*     | 145                                  | -312.99                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *cotY* — *cotX*     | 154                                  | -213.83                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *yjoB* — *rapA*     | 147                                  | -380.85                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+
    | *ptsI* — *splA*     | 93                                   | -291.13                                 | NOP     |
    +---------------------+--------------------------------------+-----------------------------------------+---------+

    --------------

 表`16.1 `__ 列出了枯草芽孢杆菌的一些基因对，这些基因的操纵子结构已知。 让我们根据表中的这些数据来计算其logistic回归模型：

.. code:: verbatim

    &gt;&gt;&gt; from Bio import LogisticRegression
    &gt;&gt;&gt; xs = [[-53, -200.78],
              [117, -267.14],
              [57, -163.47],
              [16, -190.30],
              [11, -220.94],
              [85, -193.94],
              [16, -182.71],
              [15, -180.41],
              [-26, -181.73],
              [58, -259.87],
              [126, -414.53],
              [191, -249.57],
              [113, -265.28],
              [145, -312.99],
              [154, -213.83],
              [147, -380.85],
              [93, -291.13]]
    &gt;&gt;&gt; ys = [1,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0]
    &gt;&gt;&gt; model = LogisticRegression.train(xs, ys)

这里，``xs``和``ys``是训练数据：``xs``包含每个基因对的预测变量，``ys``指定是否这个基因对属于相同操纵子（``1``，类别OP）或不同操纵子（``0``，类别NOP）。Logistic回归模型结果存储在model中，包含权重β\ :sub:`0`, β\ :sub:`1`, and β\ :sub:`2`:

.. code:: verbatim

    &gt;&gt;&gt; model.beta
    [8.9830290157144681, -0.035968960444850887, 0.02181395662983519]

注意β\ :sub:`1`是负的，这是因为具有更短基因间距离的基因对有更高的概率属于相同操纵子（类别OP）。另一方面，β\ :sub:`2`为正，因为属于相同操纵子的基因对通常有更高的基因表达轮廓相似性得分。参数β\ :sub:`0`是正值是因为在这个训练数据中操纵子基因对占据大多数。

函数``train``有两个可选参数：``update_fn``和``typecode``。``update_fn``可用来指定一个回调函数，以迭代数和对数似然值做参数。在这个例子中，我们可以使用这个回调函数追踪模型计算（使用Newton-Raphson迭代来最大化logistic回归模型的对数似然函数）进度：

.. code:: verbatim

    &gt;&gt;&gt; def show_progress(iteration, loglikelihood):
            print "Iteration:", iteration, "Log-likelihood function:", loglikelihood
    &gt;&gt;&gt;
    &gt;&gt;&gt; model = LogisticRegression.train(xs, ys, update_fn=show_progress)
    Iteration: 0 Log-likelihood function: -11.7835020695
    Iteration: 1 Log-likelihood function: -7.15886767672
    Iteration: 2 Log-likelihood function: -5.76877209868
    Iteration: 3 Log-likelihood function: -5.11362294338
    Iteration: 4 Log-likelihood function: -4.74870642433
    Iteration: 5 Log-likelihood function: -4.50026077146
    Iteration: 6 Log-likelihood function: -4.31127773737
    Iteration: 7 Log-likelihood function: -4.16015043396
    Iteration: 8 Log-likelihood function: -4.03561719785
    Iteration: 9 Log-likelihood function: -3.93073282192
    Iteration: 10 Log-likelihood function: -3.84087660929
    Iteration: 11 Log-likelihood function: -3.76282560605
    Iteration: 12 Log-likelihood function: -3.69425027154
    Iteration: 13 Log-likelihood function: -3.6334178602
    Iteration: 14 Log-likelihood function: -3.57900855837
    Iteration: 15 Log-likelihood function: -3.52999671386
    Iteration: 16 Log-likelihood function: -3.48557145163
    Iteration: 17 Log-likelihood function: -3.44508206139
    Iteration: 18 Log-likelihood function: -3.40799948447
    Iteration: 19 Log-likelihood function: -3.3738885624
    Iteration: 20 Log-likelihood function: -3.3423876581
    Iteration: 21 Log-likelihood function: -3.31319343769
    Iteration: 22 Log-likelihood function: -3.2860493346
    Iteration: 23 Log-likelihood function: -3.2607366863
    Iteration: 24 Log-likelihood function: -3.23706784091
    Iteration: 25 Log-likelihood function: -3.21488073614
    Iteration: 26 Log-likelihood function: -3.19403459259
    Iteration: 27 Log-likelihood function: -3.17440646052
    Iteration: 28 Log-likelihood function: -3.15588842703
    Iteration: 29 Log-likelihood function: -3.13838533947
    Iteration: 30 Log-likelihood function: -3.12181293595
    Iteration: 31 Log-likelihood function: -3.10609629966
    Iteration: 32 Log-likelihood function: -3.09116857282
    Iteration: 33 Log-likelihood function: -3.07696988017
    Iteration: 34 Log-likelihood function: -3.06344642288
    Iteration: 35 Log-likelihood function: -3.05054971191
    Iteration: 36 Log-likelihood function: -3.03823591619
    Iteration: 37 Log-likelihood function: -3.02646530573
    Iteration: 38 Log-likelihood function: -3.01520177394
    Iteration: 39 Log-likelihood function: -3.00441242601
    Iteration: 40 Log-likelihood function: -2.99406722296
    Iteration: 41 Log-likelihood function: -2.98413867259

一旦对数似然函数得分增加值小于0.01，迭代将终止。如果在500次迭代后还没有到达收敛，``train``函数返回并抛出一个``AssertionError``。

可选的关键字``typecode``几乎可以一直忽略。这个关键字允许用户选择要使用的数值矩阵类型。当为了避免大数据计算的内存问题时，可能有必要使用单精度浮点数（Float8，Float16等等）而不是默认的double型。

16.1.3 使用logistic回归模型进行分类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

调用``classify``函数可以进行分类。给定一个logistic回归模型和*x*\ :sub:`1`和*x*\ :sub:`2`的值（例如，操纵子结构未知的基因对），``classify``函数返回``1``或``0``，分别对应类别OP和NOP。例如，考虑基因对 *yxcE*, *yxcD* and *yxiB*, *yxiA*:

    --------------

    +-------------------------------------------------------------+
    | 表16.2：操纵子状态未知的相邻基因对   |
    +-------------------------------------------------------------+

    +-------------------+------------------------------------+---------------------------------------+
    | Gene pair         | Intergene distance *x*\ :sub:`1`   | Gene expression score *x*\ :sub:`2`   |
    +-------------------+------------------------------------+---------------------------------------+
    | *yxcE* — *yxcD*   | 6                                  | -173.143442352                        |
    +-------------------+------------------------------------+---------------------------------------+
    | *yxiB* — *yxiA*   | 309                                | -271.005880394                        |
    +-------------------+------------------------------------+---------------------------------------+

    --------------

Logistic回归模型预测*yxcE*, *yxcD*属于相同操纵子（类别OP），而*yxiB*,*yxiA*属于相同操纵子:


.. code:: verbatim

    &gt;&gt;&gt; print "yxcE, yxcD:", LogisticRegression.classify(model, [6,-173.143442352])
    yxcE, yxcD: 1
    &gt;&gt;&gt; print "yxiB, yxiA:", LogisticRegression.classify(model, [309, -271.005880394])
    yxiB, yxiA: 0

（这个结果和生物学文献报道的一致）。

为了确定这个预测的可信度，我们可以调用calculate函数来获得类别OP和NOP的概率(equations(`16.2 `__) and `16.3 `__)。对于*yxcE*, *yxcD*我们发现 

.. code:: verbatim

    &gt;&gt;&gt; q, p = LogisticRegression.calculate(model, [6,-173.143442352])
    &gt;&gt;&gt; print "class OP: probability =", p, "class NOP: probability =", q
    class OP: probability = 0.993242163503 class NOP: probability = 0.00675783649744

对于 *yxiB*, *yxiA*

.. code:: verbatim

    &gt;&gt;&gt; q, p = LogisticRegression.calculate(model, [309, -271.005880394])
    &gt;&gt;&gt; print "class OP: probability =", p, "class NOP: probability =", q
    class OP: probability = 0.000321211251817 class NOP: probability = 0.999678788748

为了确定回归模型的预测精确性，我们可以把模型应用到训练数据上：

.. code:: verbatim

    &gt;&gt;&gt; for i in range(len(ys)):
            print "True:", ys[i], "Predicted:", LogisticRegression.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0

表示除一个基因对外所有预测都是正确的。Leave-one-out分析可以对预测精确性给出一个更可信的估计，这是通过从训练数据中移除要预测的基因，再重新计算模型实现：

.. code:: verbatim

    &gt;&gt;&gt; for i in range(len(ys)):
            model = LogisticRegression.train(xs[:i]+xs[i+1:], ys[:i]+ys[i+1:])
            print "True:", ys[i], "Predicted:", LogisticRegression.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0

Leave-one-out分析显示这个logistic回归模型的预测只对两个基因对不正确，对应预测精确度为88%。

16.1.4 Logistic回归，线性判别分析和支持向量机
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Logistic回归模型类似于线性判别分析。在线性判别分析中，类别概率同样可由方程(`16.2 `__) and (`16.3 `__)给出。但是，不是直接估计系数β，我们首先对预测变量*x*拟合一个正态分布。然后通过这个正态分布的平均值和方差计算系数β。如果*x*的分布确实是正态的，线性判别分析将比logistic回归模型有更好的性能。另一方面，logistic回归模型对于偏态到正态的广泛分布更加强健。

另一个相似的方法是应用线性核函数的支持向量机。这样的SVM也使用一个预测变量的线性组合，但是是从靠近类别之间的边界区域的预测变量*x*来估计系数β。如果logistic回归模型(equations (`16.2 `__) and (`16.3 `__))能够很好的描述远离边界区域的*x*，我们可以期望logistic回归模型优于线性核函数SVM，因为它应用了更多数据。如果不是，SVM可能更好。

Trevor Hastie, Robert Tibshirani, and Jerome Friedman: *The Elements of
Statistical Learning. Data Mining, Inference, and Prediction*. Springer
Series in Statistics, 2001. Chapter 4.4.

16.2 *k*-最近邻居（KNN）
---------------------------

16.2.1 背景和目的
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*k*-最近邻居方法是一种不需要对数据训练模型的监督学习方法。数据点是基于训练数据集的*k*个最近邻居类别进行分类的。

在Biopython中，KNN方法可在``Bio.KNN``中获得。我们使用同样的操纵子数据集(`16.1 `__)来说明Biopython中KNN方法的用法。

16.2.2 初始化一个*K*NN模型
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

使用表`16.1 `__中的数据，我们创建和初始化一个*K*NN模型

.. code:: verbatim

    &gt;&gt;&gt; from Bio import kNN
    &gt;&gt;&gt; k = 3
    &gt;&gt;&gt; model = kNN.train(xs, ys, k)

这里``xs``和``ys``和`16.1.2 `__中的相同。``k``是分类中的邻居数*k*。对于二分类，为*k*选择一个奇数可以避免tied votes。函数名``train``在这里有点不合适，因为就没有训练模型：这个函数仅仅是用来存储变量``xs``，``ys``和``k``。

16.2.3 使用KNN模型来分类
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

应用KNN模型对新数据进行分类，我们使用``classify``函数。这个函数以一个数据点(*x*\ :sub:`1`,\ *x*\ :sub:`2`)）为参数并在训练数据集xs中寻找*k*-最近邻居。然后基于在这*k*个邻居中出现最多的类别（``ys``）来对数据点(*x*\ :sub:`1`,\ *x*\ :sub:`2`)进行分类。

对于基因对*yxcE*, *yxcD* and *yxiB*, *yxiA*的例子，我们发现：

.. code:: verbatim

    &gt;&gt;&gt; x = [6, -173.143442352]
    &gt;&gt;&gt; print "yxcE, yxcD:", kNN.classify(model, x)
    yxcE, yxcD: 1
    &gt;&gt;&gt; x = [309, -271.005880394]
    &gt;&gt;&gt; print "yxiB, yxiA:", kNN.classify(model, x)
    yxiB, yxiA: 0

和logistic回归模型一致，*yxcE*, *yxcD*被归为一类（类别OP），*yxiB*, *yxiA*属于不同操纵子。

函数``classify``可以指定距离函数和权重函数作为可选参数。距离函数影响作为最近邻居的*k*个类别的选择，因为这些到查询点（*x*，*y*）有最小距离的类别被定义为邻居。默认使用欧几里德距离。另外，我们也可以如示例中的使用曼哈顿距离

.. code:: verbatim

    &gt;&gt;&gt; def cityblock(x1, x2):
    ...    assert len(x1)==2
    ...    assert len(x2)==2
    ...    distance = abs(x1[0]-x2[0]) + abs(x1[1]-x2[1])
    ...    return distance
    ...
    &gt;&gt;&gt; x = [6, -173.143442352]
    &gt;&gt;&gt; print "yxcE, yxcD:", kNN.classify(model, x, distance_fn = cityblock)
    yxcE, yxcD: 1

权重函数可以用于权重投票。例如，相比于相邻较远的邻居，我们可能想给更近的邻居一个更高的权重：

.. code:: verbatim

    &gt;&gt;&gt; def weight(x1, x2):
    ...    assert len(x1)==2
    ...    assert len(x2)==2
    ...    return exp(-abs(x1[0]-x2[0]) - abs(x1[1]-x2[1]))
    ...
    &gt;&gt;&gt; x = [6, -173.143442352]
    &gt;&gt;&gt; print "yxcE, yxcD:", kNN.classify(model, x, weight_fn = weight)
    yxcE, yxcD: 1

默认所有邻居有相同权重。

为了确定这些预测的置信度，我们可以调用函数``calculate``来计算分配到类别OP和NOP的总权重。对于默认的加权方案，这样减少了每个分类的邻居数量。对于*yxcE*, *yxcD*, 我们发现

.. code:: verbatim

    &gt;&gt;&gt; x = [6, -173.143442352]
    &gt;&gt;&gt; weight = kNN.calculate(model, x)
    &gt;&gt;&gt; print "class OP: weight =", weight[0], "class NOP: weight =", weight[1]
    class OP: weight = 0.0 class NOP: weight = 3.0

这意味着``x1``，``x2``的所有三个邻居都属于NOP类别。对另一个例子*yesK*, *yesL*我们发现

.. code:: verbatim

    &gt;&gt;&gt; x = [117, -267.14]
    &gt;&gt;&gt; weight = kNN.calculate(model, x)
    &gt;&gt;&gt; print "class OP: weight =", weight[0], "class NOP: weight =", weight[1]
    class OP: weight = 2.0 class NOP: weight = 1.0

这意思是两个邻居是操纵子对，另一个是非操纵子对。

对于KNN方法的预测精确性，我们对训练数据应用模型：

.. code:: verbatim

    &gt;&gt;&gt; for i in range(len(ys)):
            print "True:", ys[i], "Predicted:", kNN.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0

显示除了两个基因对所有预测都是正确的。Leave-one-out分析可以对预测精确性给出一个更可信的估计，这是通过从训练数据中移除要预测的基因，再重新计算模型实现：

.. code:: verbatim

    &gt;&gt;&gt; for i in range(len(ys)):
            model = kNN.train(xs[:i]+xs[i+1:], ys[:i]+ys[i+1:])
            print "True:", ys[i], "Predicted:", kNN.classify(model, xs[i])
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 1
    True: 1 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 1
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 0
    True: 0 Predicted: 1

Leave-one-out分析显示这个KNN模型的预测正确17个基因对中的13个，对应预测精确度为76%。

16.3 Naive贝叶斯
-----------------

这部分将描述模块``Bio.NaiveBayes``.

16.4 最大熵
---------------------

这部分将描述模块``Bio.MaximumEntropy``.

16.5 马尔科夫模型
-------------------

这部分将描述模块``Bio.MarKovModel``和``Bio.HMM.MarKovModel``.

</article>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2013 <span title="0.33068s from fe3.rs.github.com">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/julyfire/Biopython-cn/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-remove-close close ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>

    
  </body>
</html>

