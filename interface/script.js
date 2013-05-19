/**
 * Hide variable's value
 */
function hide(elm) {
  name = elm
  elm = elm.parentNode.nextSibling;

  while(!(elm instanceof HTMLElement)) {
    elm = elm.nextSibling;
  } 

  if(elm.style.display == 'inline' || !elm.style.display) {
    elm.style.display = 'none';
    name.style.color = 'red';
  }
  else {
    elm.style.display = 'inline';
    name.style.color = 'blue';
  }
}
